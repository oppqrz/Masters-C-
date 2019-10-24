// INITIAL ATTEMPT A MINICING THE MASTER CODE IN C++ FOR IMPROVED SPEED
// THIS ATTEMPT USES ONLY 1 MASS OPTION 1 SOLAR MASS



#include <math.h>
#include <fstream>
#include <vector>
#include "Math/Polynomial.h"
#include "Math/Interpolator.h"


// GENERAL CONSTANTS
#define G 6.67e-8
#define PI 3.141562

#define Myr_seconds  1e6 * 365.25 *24*3600
#define yr_seconds  365.25 *24*3600

//SOLAR CONSTANTS
#define moi_s 7e53
#define W_s 2.8e-6
#define P_s 2*PI/W_s
#define t_s 4650*Myr_seconds
#define k_s 6.7e30
#define M_s 2e33
#define R_s 6.957e10
#define tau_s 12.9
#define L_s 3.839e33
#define v_esc_s sqrt(2*G*M_s/R_s)
#define rossby_s 2.01
#define rossby_sat 0.1


//SIMULATION CONSTANTS
#define accuracy 0.1
#define chi 10
#define mass 1
#define velocity_group 25
#define resolution 40
#define file_lines 424
#define W_0 3.0 * W_s
#define T_0 4.0 * Myr_seconds
#define T_max t_s

using namespace std;

double chi_calculation(double observation,double model){
  return pow((log10(observation)- log10(model)),2);
};

double turnover(double temp){
  return 0.002 +314.24*exp(-temp/1952.5 - pow((temp/6250),18));
};

double rossby(double w, double temp){
  double answer = ((2*PI)/(w*3600*24))/turnover(temp);
  //cout << "Rossby Number: " << answer << endl;
  return answer;
};

double mdot(double t, double lum, double rad, double temp){
  const double C = 8.5e-14;
  const double mass_SI = (mass*M_s);

  double massloss = C*lum*rad/mass_SI;
  massloss *= (1+(G*mass_SI/pow(rad,2))/4300*9.81);

  massloss *= pow(temp/4000,3.5);

  massloss *= M_s/(365.25*24*3600);

  return massloss;
};

double openflux(double rossby,double openflux_sat, double power_factor){
  //cout << rossby_sat << endl;
  if(rossby <= rossby_sat)
  {
    //cout <<"Saturated" << endl;
    return openflux_sat;
  }
  else
  {
    //cout << "Unsaturated: " << openflux_sat*pow((rossby/rossby_sat),-power_factor) << endl;
    return openflux_sat*pow((rossby/rossby_sat),-power_factor);
  };
};

double opentorque(double w, double grad,double dt, double t,double rad, double temp, double lum, double ross, double OFS, double alpha){
  const double mass_SI = (mass*M_s);
  const double k3 = 0.65;
  const double k4 = 0.06;
  const double m = 0.31;
  double ang = (w + grad*dt);

  double f = (ang*pow(rad,1.5))/sqrt(G*mass_SI);
  double v_esc = sqrt(2*G*mass_SI/rad);


  double mass_loss = mdot(t,lum,rad,temp);
  //cout << "Mass Loss: " << mass_loss << endl;
  double open_flux = openflux(ross,OFS,alpha);
  double exponent1 = (1-2*m);
  double exponent2 = (2-4*m);
  double exponent3 = 2*m;

  double frac1 = pow(open_flux,2);
  double frac2 = v_esc*sqrt(1+pow(f/k4,2));

  double torque = -pow(mass_loss,exponent1);
  torque = torque * ang * pow(rad,exponent2) * pow(k3,2);
  torque *= pow(frac1/frac2, exponent3);
  //cout <<"Torque: "<< torque <<endl;
  return torque;
};

void rk(vector<double>* w, vector<double>* t,double* T,double* moi, double* R, double* teff, double* L, double OFS, double alpha){
  ROOT::Math::Interpolator Inertia(file_lines);
  Inertia.SetData(file_lines,T,moi);

  ROOT::Math::Interpolator Radius(file_lines);
  Radius.SetData(file_lines,T,R);

  ROOT::Math::Interpolator Tempeture(file_lines);
  Tempeture.SetData(file_lines,T,teff);

  ROOT::Math::Interpolator Luminosity(file_lines);
  Luminosity.SetData(file_lines,T,L);


  int i = 0;

  while(t->at(i) < T_max)
  {
    double MoI = Inertia.Eval(t->at(i));
    //cout <<"MOI: " << MoI << endl;
    double dMoI = Inertia.Deriv(t->at(i));
    //cout << "dMOI: " <<dMoI << endl;
    double rad = Radius.Eval(t->at(i));
    //cout <<"RADIUS: " << rad << endl;
    double temp = Tempeture.Eval(t->at(i));
    //cout <<"TEMP: "<< temp << endl;
    double lum = Luminosity.Eval(t->at(i));
    //cout <<"LUM: " << lum << endl;
    double ros = rossby(w->at(i),temp);
    //cout <<"ROSSBY: "<< ros << endl;

    double torque0 = opentorque(w->at(i),0,0,t->at(i),rad,temp,lum,ros,OFS,alpha);
    double dt = abs(accuracy*MoI/(torque0/w->at(i) - dMoI));
    if(dt > 1e15){
      dt = 1e15;
    }

    double torque1 = opentorque(w->at(i),0,dt,t->at(i),rad,temp,lum,ros,OFS,alpha);
    double k1 = torque1/MoI - (w->at(i)/MoI)*dMoI;

    double torque2 = opentorque(w->at(i),k1,dt/2,t->at(i),rad,temp,lum,ros,OFS,alpha);
    double k2 = torque2/MoI - (w->at(i)/MoI)*dMoI;

    double torque3 = opentorque(w->at(i),k2,dt/2,t->at(i),rad,temp,lum,ros,OFS,alpha);
    double k3 = torque3/MoI - (w->at(i)/MoI)*dMoI;

    double torque4 = opentorque(w->at(i),k3,dt,t->at(i),rad,temp,lum,ros,OFS,alpha);
    double k4 = torque4/MoI - (w->at(i)/MoI)*dMoI;

    double kfinal = (k1+2*k2+2*k3+k4)/6;

    w->push_back(w->at(i) + (kfinal*dt));
    t->push_back(t->at(i) + dt);

    //cout << " Age: " << t->at(i) + dt;
    //cout << " Size: " << w.size() << endl;
    i +=1;
  };

return;
};




void FasterMaster(){
  cout << "Start" << rossby_sat << endl;

  time_t timer;
  //TCanvas* canvas = new TCanvas("canvas")
  double startTime = time(&timer);

  ifstream structure_file;

  double M[file_lines];
  double T[file_lines];
  double teff[file_lines];
  double L[file_lines];
  double R[file_lines];
  double k2con[file_lines];
  double k2rad[file_lines];
  double moi[file_lines];

  float skip;

  structure_file.open("1SolarMass.txt");
  for(int i = 0; i < file_lines; i++)
  {
    structure_file >> M[i];
    structure_file >> T[i];
    structure_file >> teff[i];
    structure_file >> L[i];
    structure_file >> skip;
    structure_file >> R[i];
    structure_file >> skip;
    structure_file >> skip;
    structure_file >> skip;
    structure_file >> skip;
    structure_file >> skip;
    structure_file >> k2con[i];
    structure_file >> k2rad[i];

    M[i] = M[i] * M_s;
    T[i] = pow(10,T[i])*yr_seconds;
    L[i] = pow(10,L[i]) * L_s;
    R[i] = R[i] * R_s;
    moi[i] = (k2con[i]*k2con[i]+k2rad[i]*k2rad[i])*M[i]*R[i]*R[i];


  };
  structure_file.close();

  const int data_lines = 14;
  double data_t[data_lines];
  double slow_v[data_lines];
  double slow_er[data_lines];
  double medium_v[data_lines];
  double medium_er[data_lines];
  double fast_v[data_lines];
  double fast_er[data_lines];

  ifstream data_file;
  data_file.open("Gbdata_chi_M1.txt");
  for(int i = 0; i < data_lines; i++)
  {
    data_file >> data_t[i];
    data_file >> slow_v[i];
    data_file >> slow_er[i];
    data_file >> medium_v[i];
    data_file >> medium_er[i];
    data_file >> fast_v[i];
    data_file >> fast_er[i];
  };
  data_file.close();

  double x_start = 1200e18;//5e19;
  double x_finish = 5000e18;//1e21;
  double x_inc = (x_finish - x_start)/resolution;
  double y_start = 1.8;
  double y_finish = 2.8;
  double y_inc = (y_finish - y_start)/resolution;

  vector<double> x_range;
  for(int i = 0; i <= resolution; i++)
  {
    x_range.push_back(x_start + x_inc*i);
    //cout << x_range[i] << endl;
  }

  vector<double> y_range;
  for(int i = 0; i <= resolution; i++)
  {
    y_range.push_back(y_start + y_inc*i);
    //cout << y_range[i] << endl;
  }

  //  This section of the code greats heatmaps

  TCanvas* canvas = new TCanvas("Canvas","Parameter Fitting",1000,1000);
  //canvas->SetLogx();
  //canvas->SetLogy();
  TH2D* histogram = new TH2D("histogram","Chi Squared",resolution+1,x_start-1e-6,x_finish+1e6,resolution+1,y_start-1e-6,y_finish+1e-6);
  double chisquared = 0;

  cout << "Doing stuff" << endl;

  for(int i_OFS = 0; i_OFS <= resolution ; i_OFS++)
  {
      //cout << "Outer " << i_OFS << endl;

      for(int i_alpha = 0 ; i_alpha <= resolution ; i_alpha++)
      {
        //cout << i_alpha << endl;

        vector<double>* w = new vector<double>() ;
        vector<double>* t = new vector<double>() ;
        w->push_back(W_0);
        t->push_back(T_0);

        double OFS = x_range[i_OFS];
        double alpha = y_range[i_alpha];
        //cout << "X: " << OFS << " Y: " << alpha << endl;

        rk(w,t,T,moi,R,teff,L,OFS,alpha);

        for(int i = 0 ;i < t->size();i++)
        {
          t->at(i) /= Myr_seconds;
          w->at(i) /= W_s;
        };

        double* sim_t = &t->at(0);
        double* sim_w = &w->at(0);

        ROOT::Math::Interpolator simulation(t->size());
        simulation.SetData(t->size(),sim_t,sim_w);

        chisquared = 0;

        for(int j =0; j < data_lines ; j++)
        {
            if(data_t[j] < T_0/Myr_seconds)
          {
            continue;
          } else
          {
            chisquared += chi_calculation(slow_v[j],simulation.Eval(data_t[j]));
          }
        };

        chisquared /= double(data_lines);
        //cout << "x: "<< OFS << " y: " << alpha << " z: " << 1/chisquared << endl;

        histogram->Fill(OFS,alpha,1/chisquared);
        //histogram->Fill(OFS,alpha);

      };
      cout << i_OFS <<"%" << endl;

  };

  histogram->Draw("surf2");

  double endTime = time(&timer);
  double runTime = endTime-startTime;
  cout << "RUN TIME:" << runTime << endl;


  double endTime = time(&timer);
  double runTime = endTime-startTime;
  cout << "RUN TIME:" << runTime << endl;

};
