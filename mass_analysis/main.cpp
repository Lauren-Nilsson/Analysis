#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <boost/lexical_cast.hpp>
#include "oligomer.h"
#include "vector3d.h"
#include "bead.h"
#include "subunit.h"

using namespace std;


namespace patch
{
   template < typename T > std::string to_string( const T& n )
   {
      std::ostringstream stm ;
      stm << n ;
      return stm.str() ;
   }
}

class RunningStat //Welford's algorithm from https://www.johndcook.com/blog/standard_deviation/
{
public:
   RunningStat() : m_n(0) {}
   
   void Clear()
   {
      m_n = 0;
   }
   
   void Push(double x)
   {
      m_n++;
      
      // See Knuth TAOCP vol 2, 3rd edition, page 232
      if (m_n == 1)
      {
         m_oldM = m_newM = x;
         m_oldS = 0.0;
      }
      else
      {
         m_newM = m_oldM + (x - m_oldM)/m_n;
         m_newS = m_oldS + (x - m_oldM)*(x - m_newM);
         
         // set up for next iteration
         m_oldM = m_newM; 
         m_oldS = m_newS;
      }
   }
   
   int NumDataValues() const
   {
      return m_n;
   }
   
   double Mean() const
   {
      return (m_n > 0) ? m_newM : 0.0;
   }
   
   double Variance() const
   {
      return ( (m_n > 1) ? m_newS/(m_n - 1) : 0.0 );
   }
   
   double StandardDeviation() const
   {
      return sqrt( Variance() );
   }
   
private:
   int m_n;
   double m_oldM, m_newM, m_oldS, m_newS;
};



VECTOR3D dist(BEAD* A, BEAD* B){                            //finds distance considering periodic boundaries.
    VECTOR3D r_vec; //= (A->pos - B->pos);
    r_vec.x = A->pos.x - B->pos.x;
    r_vec.y = A->pos.y - B->pos.y;
    r_vec.z = A->pos.z - B->pos.z;
    VECTOR3D box = A->bx;
    if (r_vec.x>box.x/2) r_vec.x -= box.x;
    if (r_vec.x<-box.x/2) r_vec.x += box.x;
    if (r_vec.y>box.y/2) r_vec.y -= box.y;
    if (r_vec.y<-box.y/2) r_vec.y += box.y;
    if (r_vec.z>box.z/2) r_vec.z -= box.z;
    if (r_vec.z<-box.z/2) r_vec.z += box.z;
    return r_vec;
}


double compute_MD_trust_factor_R(int &hiteqm, bool &done, string directory) {
  
   string inPath = directory+"/energy.out";
   ifstream in(inPath.c_str(), ios::in);
   if (!in) {
      if (!in)
         cout << "ERR: FILE " << directory+"/energy.out" << " NOT OPENED. Check directory and/or filename." << endl;
      return 0;
   }
   string dummy;
   string col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11;
   double col2_db, col7_db, col8_db, col5_db, col1_db;
   vector<double> ext, ke, pe, lj, time;
   in >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
   while (in >> col1 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7 >> col8 >> col9) {
      try {  
         col1_db = boost::lexical_cast<double>(col1);
         col2_db = boost::lexical_cast<double>(col2); //use boost lexical_cast to look for extra headers from restart files
         col7_db = boost::lexical_cast<double>(col7);
         col8_db = boost::lexical_cast<double>(col8);
         col5_db = boost::lexical_cast<double>(col5);
      } catch (boost::bad_lexical_cast&) {
         cout << "Caught a restart header!" << endl;
         goto next;
      }
      in >> col10 >> col11;
      time.push_back(col1_db);
      ext.push_back(col7_db);
      ke.push_back(col2_db);
      pe.push_back(col8_db);
      lj.push_back(col5_db);
      next: ;
      
   }
   //determine when equilibrium is reached
   int energyIndex = 0;
   double x_sum, y_sum, xy_sum, x_m, y_m;
   double xx_sum, slope;
   int number_sections = 10;
   int section_size = floor(pe.size()/number_sections);
   vector<double> slope_vec;
   slope_vec.resize(number_sections);
   vector<double> stdev_vec;
   stdev_vec.resize(number_sections);
   cout << "Each section has " << section_size << " points." << endl;
   
   vector<RunningStat> RS;
   RS.resize(number_sections);
   double mean;
   double variance;
   double stdev;
   
   for (int i = 0; i < number_sections; i++){ //loop over sections
      x_sum = 0;
      y_sum = 0;
      xy_sum = 0;
      xx_sum = 0;
      for (int j = 0; j < section_size; j++){ //loop over points in sections
         x_sum += j;
         y_sum += pe[(i*section_size + j)];
         xy_sum += j * pe[(i*section_size + j)];
         xx_sum += j * j;
         
         RS[i].Push(pe[(i*section_size + j)]);
      }
      x_m = x_sum / double(section_size);
      y_m = y_sum / double(section_size);
      slope = ((xy_sum) - (double(section_size) * x_m * y_m) ) / (xx_sum - (double(section_size) * x_m * y_m) ); 
   //   cout << "For section " << i << " the slope is " << slope << endl;
      slope_vec[i] = slope;
      
      mean = RS[i].Mean();
      variance = RS[i].Variance();
      stdev = RS[i].StandardDeviation();
     // cout << "standard deviation is " << stdev << endl << endl;
      
      stdev_vec[i] = stdev;
      
   }
   
   RunningStat rs;
   
   
   
   
   //find sections which have a slope +- 1e-6
   vector<bool> flat_vec;
   flat_vec.resize(number_sections);
   for (int i = 0; i < number_sections; i++) {
      //if (slope_vec[i] < 1e-6 && slope_vec[i] > -1e-6) {
      if (stdev_vec[i] < 0.02 && slope_vec[i] < 1e-4 && slope_vec[i] > -1e-4) {
         flat_vec[i] = true;
       //  cout << "section " << i << " is flat." << endl; 
      }
      else flat_vec[i] = false;
   }
   flat_vec[9] = true;
   //See how many of the final sections are at equilibrium
   int flat_sections = 0;
   for (int i = (number_sections - 1); i > -1; i--) {
      if (flat_vec[i] == true) {
         flat_sections += 1;
         hiteqm = i*section_size;
      }
      else break;
   }
   hiteqm = time[hiteqm];
   cout << "There are " << flat_sections << " ending consecutive flat sections" << endl;
   done = true;
   //Let the user know how long equilibrium has been reached (or if it needs restart)
   if (flat_sections > 0) {
      cout << 100.0*double(flat_sections)/double(number_sections) << "% of simulation is in equilibrium. Analyzing..." << endl;
      cout << "Equilibrium starts at timestep " << hiteqm << "." << endl;
   } else if (flat_sections ==0) {
      cout << "Simulation has not reached equilibrium! Analysis aborted." << endl;
      done = false;
   } 
   
   //Compute R
   double ext_mean = 0;
   for (unsigned int i = 0; i < ext.size(); i++)
      ext_mean += ext[i];
   ext_mean = ext_mean / (ext.size() - 0 );
   double ke_mean = 0;
   for (unsigned int i = 0; i < ke.size(); i++)
      ke_mean += ke[i];
   ke_mean = ke_mean / ke.size();
   
   double ext_sd = 0;
   for (unsigned int i = 0; i < ext.size(); i++)
      ext_sd += (ext[i] - ext_mean) * (ext[i] - ext_mean);
   ext_sd = ext_sd / (ext.size() - 0 );
   ext_sd = sqrt(ext_sd);
   
   double ke_sd = 0;
   for (unsigned int i = 0; i < ke.size(); i++)
      ke_sd += (ke[i] - ke_mean) * (ke[i] - ke_mean);
   ke_sd = ke_sd / (ke.size() - 0 );
   ke_sd = sqrt(ke_sd);
   
   double R = ext_sd / ke_sd;
   //    if (world.rank() == 0)
   //    {
   ofstream out( (directory+"/analysis.rdat").c_str() );
   out << "Sample size " << ext.size() << endl;
   out << "Sd: ext, kinetic energy and R" << endl;
   out << ext_sd << setw(15) << ke_sd << setw(15) << R << endl;
   //    }
   cout << endl << endl << "R is: " << R << endl << endl;
   
   return R;
}



int main() {
   
   { std::cout << __cplusplus << std::endl; }
   
   //cout << "Please enter the directory you want to analyze." << endl;
   string old_directory;
   //getline (cin, directory);
   //old_directory = "/run/media/lm44/Extreme SSD/DATA/CG-T1/300_R2_0.002/";
   //old_directory = "/run/media/lm44/Extreme SSD/DATA/T3_T4/200_patterns_R1/";
   old_directory = "/run/media/lm44/Extreme SSD/DATA/T1/150_elasticity/";

   int number_of_particles;



   string directory;
   
   
   int protein_vec[] = {1000};
   int salt_vec[] = {1000};
 
   for (int p_i = 0; p_i < (sizeof(protein_vec)/sizeof(protein_vec[0])); p_i ++) {
      for (int s_i = 0; s_i < (sizeof(salt_vec)/sizeof(salt_vec[0])); s_i ++) {
         //string file_name = string("R1CGT1_sweep_41nvt_300_")+patch::to_string(protein_vec[p_i])+"_"+patch::to_string(salt_vec[s_i])+"/outfiles";
         //string file_name = string("T3_pattern_sweep_41nvt_200_")+patch::to_string(protein_vec[p_i])+"_"+patch::to_string(salt_vec[s_i])+"/outfiles";
         string file_name = string("T1_sweep_43nvt_150_300_500_")+patch::to_string(protein_vec[p_i])+"_"+patch::to_string(salt_vec[s_i])+"/outfiles";
         directory = old_directory+file_name;
         string coords = directory+"/ovito.lammpstrj";
         string nrg = directory+"/energy.out";
         string analysis = directory+"/analysis";
         
         ifstream crds;                                           //open coordinates file
         crds.open(coords.c_str());
         if (!crds) {                                             //check to make sure file is there
            cerr << "ERR: FILE " << coords << " NOT OPENED. Check directory and/or filename.";
            exit(1);
         }
         
         cout << "Mass spectrum file will be saved to " << directory+"/analysis.ms" << endl;
         cout << "Pair-correlation file will be saved to " << directory+"/analysis.gr" << endl;
         cout << "R value will be saved to " << directory+"/analysis.rdat" << endl;
         cout << "Sub-oligomer structural information will be saved to " << directory+"/analysis.so" << endl;
         
         string dummy ;                                           //dummy string for text in file
         int number_of_beads;                                     //how many beads in the file
         int number_of_subunits = 0;
         //int number_of_timesteps = 15920;                       //how many timesteps in the file
         int bead_index;                                          //index of bead
         int type;                                                //type of bead
         long double x, y, z;                                     // x y z coordinates
         long double box_size;
         
         vector<OLIGOMER> oligomers_list;                         //create vector to hold oligomers for mass spectrum analysis
         vector<OLIGOMER> sub_oligomers_list;                     //create vector to hold sub_oligomers for morphology analysis
         vector<BEAD> subunit_bead;                               //Create particles, named subunit_bead
         vector<SUBUNIT> protein;                                 //create subunits, named protein
         
         ofstream msdata( (directory+"/analysis.ms").c_str() );
         ofstream sodata( (directory+"/analysis.so").c_str() );
         ofstream grdata( (directory+"/analysis.gr").c_str() );
         ofstream rdata( (directory+"/analysis.rdat").c_str() );
         ofstream grinfo( (directory+"/info.gr").c_str() );
         
         int mstime = -1;                                         //parameter for ms_bin filling      
         vector<int> massbins(protein.size());
         vector<vector<int> > ms_bin;
         int sotime = -1;                                         //parameter for so_bin filling      
         vector<int> sobins(protein.size()*6);
         vector<vector<int> > so_bin;
         vector<int> pair_bins;
         
         int hiteqm = 0;
         bool done;
         int a = 0; //index
         double delta_g;
         int gr_size = 100; //bin number for g(r)
         int beadType = 5; // which type of bead to analyze with pcf?
         int num_beadType = 0; //number of beads in the box of beadType
         double rho; //density of beadType in the box
         
         compute_MD_trust_factor_R( hiteqm, done, directory );
         
         if (done){
            while ( crds >> dummy ){ //while it hasn't reached the end of the line yet
               int timestep;
               crds >> dummy >> timestep >> dummy >> dummy >> dummy >> dummy ;      //reading the header of the timestep
               crds >> number_of_beads;
               crds >> dummy >> dummy >> dummy >> dummy >> box_size >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy ;
               if (a==0){
                  if (number_of_beads % 41 == 0) number_of_particles = 41;
                  if (number_of_beads % 43 == 0) number_of_particles = 43;
                  subunit_bead.resize(number_of_beads);
                  number_of_subunits = number_of_beads / number_of_particles;
                  protein.resize(number_of_subunits);
                  //ms_bin.resize(number_of_timesteps, vector<int>(number_of_subunits) );   
                  for (int b = 0; b < number_of_beads; b++){
                     subunit_bead[b].bx.x = box_size*2;
                     subunit_bead[b].bx.y = box_size*2;
                     subunit_bead[b].bx.z = box_size*2;
                  }
                  delta_g = subunit_bead[0].bx.x / (2*gr_size);
               }

               int myindex;
               pair_bins.resize(gr_size);
               
               for (int b = 0; b < number_of_beads; b++){
                  crds >> bead_index >> type >> x >> y >> z >> dummy >> dummy;
                  myindex = b / (int)number_of_particles;
                  subunit_bead[b].id = bead_index-1;
                  
                  subunit_bead[b].type = type;
                  subunit_bead[b].pos.x = x;
                  subunit_bead[b].pos.y = y;
                  subunit_bead[b].pos.z = z;
                  if (a==0){
                     protein[myindex].itsB.push_back(&subunit_bead[b]);
                     protein[myindex].id = myindex;
                     if (subunit_bead[b].type == beadType) num_beadType += 1; //count how many beads are of beadType
                  }
               }
               
               if (a==0) { //determine rho
                  rho = num_beadType / (subunit_bead[0].bx.x * subunit_bead[0].bx.x * subunit_bead[0].bx.x); //number density ******** (ASSUMPTION sigma = 1nm!) *******
               }
               //  cout << "Read section: " << a+1 << " , " ;
         
         
               if (timestep >= hiteqm && (a % 5 == 0)) {
                  cout << "Analyzing position file " << a << " with timestep " << timestep << endl;
                  ms_bin.push_back(vector<int>() );
                  so_bin.push_back(vector<int>() );
                  //Generate pair_correlation bin
                  for (unsigned int i = 0; i < number_of_beads; i++){
                     for (unsigned int j = i + 1; j < number_of_beads; j++) {
                        if (subunit_bead[i].type == beadType && subunit_bead[j].type == beadType){
                           double distance = dist(&subunit_bead[i], &subunit_bead[j]).GetMagnitude();
                           if (distance < (box_size)) {
                              int gr_index = floor(distance/delta_g);
                              pair_bins[gr_index] += 1;
                           }
                        }
                     }
                  }
               // Save pair data to a file
               for (unsigned int j = 0; j < pair_bins.size(); j++) {
                  grdata << pair_bins[j] << setw(15);                         //print gr bin data to file
                  pair_bins[j] = 0;                                           //Clear gr bin
               }
               grdata << endl;
            
            
               int index = -1;
               for (unsigned int i = 0; i < protein.size(); i++) {            //Create oligomers for mass spectrum analysis
                  int oldsize = 0;
                  
                  if (protein[i].itsO.size() == 0) {                          //if the unit isn't already counted...
                     
                     oligomers_list.push_back(OLIGOMER(VECTOR3D(0, 0, 0)));   //create an oligomer for the unit
                     index += 1;
                     oligomers_list[index].itsS.push_back(&protein[i]);       //add unit to oligomer
                     oligomers_list[index].id = index;
                     protein[i].itsO.push_back(oligomers_list[index]);        //add oligomer to unit
                     while (oldsize < oligomers_list[index].itsS.size()) {    //while the oligomer is still growing...
                        int n = oldsize;
                        oldsize = (int) oligomers_list[index].itsS.size();    //see how much the oligomer has grown
                        for (int j = n; j < oldsize; j++) {                   //loop over the growth from last round
                           int g = oligomers_list[index].itsS[j]->id;
                           for (unsigned int k = i + 1; k < protein.size(); k++) { //look for new growth
                              if (protein[k].itsO.size() == 0) {              //if it isn't in an oligomer yet...
                                 for (unsigned int m = 0;
                                    m < protein[g].itsB.size(); m++) {        //check to see if it is in this oligomer
                                    for (unsigned int n = 0; n < protein[k].itsB.size(); n++) {
                                       if (dist(protein[g].itsB[m], protein[k].itsB[n]).GetMagnitude() < 3) {
                                           //  if (protein[k].itsB[n]->id == 0 && protein[g].itsB[m]->id == 0){
                                                oligomers_list[index].itsS.push_back(&protein[k]);    //if it is attached, add it
                                                protein[k].itsO.push_back(oligomers_list[index]);     //mark subunit as bonded
                                                goto finish;
                                         //    }
                                          } //if
                                       } //for n
                                    } //for m
                                 } //if
                                 finish:;
                              } //for k
                           } //for j
                        } //while
                     } //if
                  } //for i
                     
                  mstime += 1;
                  ms_bin[mstime].resize(protein.size());
                  for (unsigned int i = 0; i < oligomers_list.size(); i++) {
                     if (oligomers_list[i].itsS.size() >= 1) {
                        ms_bin[mstime][(oligomers_list[i].itsS.size() - 1)] += 1;//fill mass bins
                     }
                  }
                  
                  for (unsigned int j = 0; j < ms_bin[mstime].size(); j++) {
                     msdata << ms_bin[mstime][j] << setw(15);                    //print mass bin data to file
                  }
                  msdata << endl;
                  
                  
                  
                  
                  index = -1;
                  for (unsigned int i = 0; i < oligomers_list.size(); i++) {              //Loop over oligomers
                    // cout << "Looking at oligomer " << i << " of size " << oligomers_list[i].itsS.size() << endl;
                     for (unsigned int j = 0; j < oligomers_list[i].itsS.size(); j++){    //Loop over subunits within the oligomer
                  //      cout << "Looking at subunit " << j << endl;
                        for (unsigned int jj = 0; jj < oligomers_list[i].itsS[j]->itsB.size(); jj++){ //Loop over beads within the subunit
                        //   cout << "Looking at bead " << jj << " of subunit " << j << endl;
                           int oldsize = 0;
                           if (oligomers_list[i].itsS[j]->itsB[jj]->itsSO.size() == 0 && oligomers_list[i].itsS[j]->itsB[jj]->type == 0){ // if bead not counted yet
                           //   cout << "Looking at bead " << jj << " of subunit " << j << endl;
                              sub_oligomers_list.push_back(OLIGOMER(VECTOR3D(0,0,0)));    //create a sub_oligomer for the bead
                              index += 1;
                              sub_oligomers_list[index].itsB.push_back(oligomers_list[i].itsS[j]->itsB[jj]); //add bead to oligomer
                              sub_oligomers_list[index].id = index;
                              oligomers_list[i].itsS[j]->itsB[jj]->itsSO.push_back(sub_oligomers_list[index]); //add oligomer to bead
                          //    cout << "New oligomer created!" << endl;
                              while (oldsize < sub_oligomers_list[index].itsB.size()) {   //while the oligomer is still growing
                                 int n = oldsize;
                                 oldsize = (int) sub_oligomers_list[index].itsB.size();   //see how much the oligomer has grown
                               //  cout << "Oligomer has grown by " << oldsize - n << endl;
                                 for (int growth = n; growth < oldsize; growth++) {                      //loop over the growth from the last round
                                    int g = sub_oligomers_list[index].itsB[growth]->id;
                                  //  cout << "Looking at new bead, " << g << endl;
                                    for (unsigned int k = 0; k < oligomers_list[i].itsS.size(); k++) {     //Look for new growth
                                       for (unsigned int kk = 0; kk < oligomers_list[i].itsS[k]->itsB.size(); kk++) {   //Loop over other beads
                                          if (oligomers_list[i].itsS[k]->itsB[kk]->itsSO.size() == 0 && oligomers_list[i].itsS[k]->itsB[kk]->type == 0) {//if it isn't in an oligomer yet
                                             int newg = (int) oligomers_list[i].itsS[k]->itsB[kk]->id;
                                          //   cout << "Evaluating bead, " << newg << " for attachment" << endl;
                                             if (dist(&subunit_bead[g], &subunit_bead[newg]).GetMagnitude() < 2) { //check to see if it is in this oligomer
                                               // cout << "Found new bead! Adding " << newg << " from " << g << " to " << index << " oligomer" << endl;
                                                sub_oligomers_list[index].itsB.push_back(oligomers_list[i].itsS[k]->itsB[kk]); //if it is attached, add it
                                                oligomers_list[i].itsS[k]->itsB[kk]->itsSO.push_back(sub_oligomers_list[index]); //mark bead as bonded
                                                goto finish2;
                                             }
                                          }
                                          finish2:;
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
                   //cout << "Finished so analysis" << endl; 
            //       cout << "There are " << sub_oligomers_list.size() << " sub_oligomers" << endl;
                   
                  sotime += 1;
                 // cout << "Got here" << endl;
                  so_bin[sotime].resize(protein.size()*2);
               //   cout << "Got here!" << endl;
                  for (unsigned int i = 0; i < sub_oligomers_list.size(); i++) {
                     if (sub_oligomers_list[i].itsB.size() >= 1) {
                        so_bin[sotime][((sub_oligomers_list[i].itsB.size() - 1)) / 3] += 1;//fill so bins
                     //   cout << "Filling " << sotime << " x " << (sub_oligomers_list[i].itsB.size() - 1) << endl;
                     }
                  }
               //   cout << "Filled bins" << endl;
                  
                  for (unsigned int j = 0; j < so_bin[sotime].size(); j++) {
                     sodata << so_bin[sotime][j] << setw(15);                    //print mass bin data to file
                  }
                  sodata << endl;
                  
                  cout << "Printed to file" << endl;
 
 
                  
                  for (unsigned int i = 0; i < protein.size(); i++) {            // clear oligomer pointers from subunit
                     protein[i].itsO.clear();
                     for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
                        protein[i].itsB[ii]->itsSO.clear();
                     }
                  }
                  
                  oligomers_list.erase(oligomers_list.begin(),oligomers_list.end());               //erases oligomer objects
                  sub_oligomers_list.erase(sub_oligomers_list.begin(),sub_oligomers_list.end());               //erases oligomer objects
                  
                  
                  
                  
               } //if at equilibrium
               a++;
            } //while file
            
            //Print some extra info for pair correlation file
            grinfo << "Number of Subunits: " << number_of_subunits << endl;
            grinfo << "PCF bin size: " << delta_g << endl;
            grinfo << "PCF bin number: " << gr_size << endl;
            grinfo << "Box size: " << box_size*2 << " " << box_size*2 << " " << box_size*2 << endl;
            grinfo << "Number density: " << rho << endl;
            grinfo << "Bead Type: " << beadType << endl;
         }//if (done)
            
         if (!done) {
            remove( (directory+"/info.gr").c_str() );
         }
      } //for salt
   } // for protein
   
   return 0 ;			
} //fxn










