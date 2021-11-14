/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <math.h> 

#include "helper_functions.h"

// using std::string;
// using std::vector;
using namespace std;
default_random_engine gen;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
	// Define normal distribution for noise
  // set the number of particles (less particles, runs quicker)
  num_particles = 150;
  //resize number of particles 
  weights.resize(num_particles);
  // Create normal Gaussian distribution for x, y and theta = yaw.
  // based on sigma values, move (x,y,theat) values
  normal_distribution<double> distance_x(x, std[0]);
  normal_distribution<double> distance_y(y, std[1]);
  normal_distribution<double> distance_theta(theta, std[2]);
    
  // Initializes particles - from the normal distributions set created in the above step
  for (int i = 0; i < num_particles; i++) {
      
    // Go to each particle and generate a x,y,theta that move along the sigma based on given uncertainty. 
    Particle particle;
    particle.x = distance_x(gen);
    particle.y = distance_y(gen);
    particle.theta = distance_theta(gen);
    particle.id = i;
    particle.weight = 1;
    
    
    particles.push_back(particle);
      
  }
    // Show only initialized; no prediction needed
  is_initialized = true;

}
 

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  // use different equation based on value of yaw rate from lesson (zero or non-zero)
  for (int i = 0; i < num_particles; i++) {
    
    double updt_x;
    double updt_y;
    double updt_theta;
    
    if (fabs(yaw_rate) < 0.0001) { //(yaw_rate == 0) 
      // Formulae for adding measurements to particles when yaw is zero
      updt_x = particles[i].x+velocity*delta_t*cos(particles[i].theta);
      updt_y = particles[i].y+velocity*delta_t*sin(particles[i].theta);
      updt_theta = particles[i].theta;
     
    } else {
      // Formulae for adding measurements to particles when yaw is nonzero
      updt_x = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta + (yaw_rate * delta_t)) - sin(particles[i].theta));
      updt_y = particles[i].y + (velocity/yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate * delta_t)));
      updt_theta = particles[i].theta +yaw_rate*delta_t;
      
    }
  // Make normal (Gaussian) distributions for adding noise
    normal_distribution<double> distance_x(updt_x, std_pos[0]);
    normal_distribution<double> distance_y(updt_y, std_pos[1]);
    normal_distribution<double> distance_theta(updt_theta, std_pos[2]);

    particles[i].x = distance_x(gen);
    particles[i].y = distance_y(gen);
    particles[i].theta = distance_theta(gen);
    
  }
}




void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  // this step is implemented in updateWeight function 
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
  //clear weights 
  weights.clear();
  
  // vehicle to map coordiante system tranformation
  
  for (unsigned int i=0; i<particles.size();i++){
    vector<LandmarkObs> trans_obs;
    vector<LandmarkObs> land_pred;
    
    for (unsigned int k=0; k<observations.size(); k++){
      LandmarkObs t_o;
      t_o.id = -1;
      t_o.x = observations[k].x * cos(particles[i].theta) - observations[k].y * sin(particles[i].theta) + particles[i].x;
      t_o.y = observations[k].x * sin(particles[i].theta) + observations[k].y * cos(particles[i].theta) + particles[i].y;

      trans_obs.push_back(t_o);
    }
    
    // land_pred is a set of landmarks in sensor range 
    vector<Map::single_landmark_s> landmarks = map_landmarks.landmark_list;
    for (unsigned int k = 0; k <landmarks.size(); k++) {
      double landmark_distance = dist(particles[i].x, particles[i].y, landmarks[k].x_f,landmarks[k].y_f );
      if (landmark_distance<sensor_range){
        LandmarkObs l_p;
        
        l_p.x = landmarks[k].x_f;
        l_p.y = landmarks[k].y_f;
        
        l_p.id = landmarks[k].id_i;
        
        land_pred.push_back(l_p);
      }
    }
    // data association step - match transformed observations to predicted landmarks 
      double temp_min_distance = 9999;
      double actual_dist = 0;
      double main_map_id = -1;
      for (unsigned int j = 0; j < trans_obs.size(); j++) {
        for (unsigned int m = 0; m<land_pred.size(); m++){
          actual_dist = dist(trans_obs[j].x,trans_obs[j].y, land_pred[m].x, land_pred[m].y);
          if (actual_dist < temp_min_distance){
            temp_min_distance = actual_dist;
            main_map_id = land_pred[m].id;
          }
        }
        trans_obs[j].id = main_map_id; 
      }

    // weights assignements based on the distace between landmarks and observations 
    double p = 1;
    double mvGd;
    
    for (unsigned int l = 0; l < land_pred.size(); l++) {
      double temp_min_distance = 9999; //reset min distance to a high number 
      int main_map_id = -1;  // deafult ID 
      double pred_x = land_pred[l].x;
      double pred_y = land_pred[l].y;
      
      for (unsigned int t = 0; t < trans_obs.size(); t++) {
        double trans_x = trans_obs[t].x;
        double trans_y = trans_obs[t].y;
        
        double actual_dist = dist(pred_x, pred_y, trans_x, trans_y);
        
        if (actual_dist< temp_min_distance){
          temp_min_distance = actual_dist;
          main_map_id = t;
        }
      }
      
      if (main_map_id != -1){
        //multivariate gaussian 
      mvGd = exp(-((pred_x-trans_obs[main_map_id].x)*(pred_x-trans_obs[main_map_id].x)/(2*std_landmark[0]*std_landmark[0]) + (pred_y-trans_obs[main_map_id].y)*(pred_y-trans_obs[main_map_id].y)/(2*std_landmark[1]*std_landmark[1]))) / (2*M_PI*std_landmark[0]*std_landmark[1]);
        
        
      p = p * mvGd;
      }
    }
    
    weights.push_back(p);
    
    particles[i].weight = p;
  }


}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  // resampled particle vector
  vector<Particle> resample_p (num_particles);
  
  // random particle index generator
  random_device rd;
  default_random_engine gen(rd());
  
  
  for (int n = 0; n < num_particles; n++) {
    discrete_distribution<int> index(weights.begin(), weights.end());
    resample_p[n] = particles[index(gen)];  
  }
  particles = resample_p;
}



void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}