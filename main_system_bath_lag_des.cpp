/*
@Purpose C++ code to compute the Lagrangian descriptors for the De Leon-Berne 
model of isomerization

@author Shibabrat Naik [2019]
*/

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <algorithm> 
// #include <sstream>

using namespace std;

const int dim = 9; // number of odes + 1 for function M

// Parameters of the system (De Leon-Berne model of isomerization)
double MASS_A = 8.0, MASS_B = 8.0, EPSILON_S = 1.0, D_X = 10.0, E = 1.510;
double Vdagger = 1.0, YW = 1.0/sqrt(2.0);
double ZETA = 0.20, LAMBDA = 1.00;
// double ZETA = 1.00, LAMBDA = 1.00;
// double ZETA = 1.00, LAMBDA = 1.50;
// double ZETA = 2.30, LAMBDA = 1.95;



/* function prototypes */
struct quadratic_params{
    double E, D_X, ZETA, LAMBDA, EPSILON_S;
};

double quadratic (double x, void *params);
double quadratic_deriv (double x, void *params);
void quadratic_fdf (double x, void *params, double *y, double *dy);

double Hderiv (double x, void *params);
double Hderiv_deriv (double x, void *params);
void Hderiv_fdf (double x, void *params, double *y, double *dy);

vector<double> get_energybounds_sosatyw(double E);
double get_pxmax_sosatyw(double E, double guess_xcoord);
double get_xbound_sosatyw(double E, double guess_xcoord);

double deleonberne_get_py(double x, double y, double px, double E);
int funcM_2dof(int xRes, int yRes, vector<double> xpx_energybounds, double tau);
double get_traj_funcM(vector<double> time_interval, vector<double> init_cond, int file_id);
double gsl_integrate(vector<double> time_interval, vector<double> init_cond, double delta_t, int file_id);
int gsl_LD_escape(vector<double> time_interval, vector<double> init_cond, double delta_t, int file_id, double& Mt, double& escape_time);
int write_matrix(const vector<vector<double> >& trajectory, double tmax, double mu, int file_id);
int deleonberne_rhs(double t, const double posIn[], double velOut[], void *params);
int deleonberne_bath_rhs(double t, const double posIn[], double velOut[], void *params);


int main(int argc, char *argv[]){

    if (argc == 1){
        puts("Computing Lagrangian descriptor for DEFAULT resolution and integration time");
        int xRes = 100;
        int yRes = 100;
        double tau = 10.0;
        puts("Calculating bounds of energy surface's intersection on the (x,px) section ...");
        
        vector<double> xpx_energybounds = get_energybounds_sosatyw(E);
        // cout << xpx_energybounds[3] << endl;

        int succ  = funcM_2dof(xRes, yRes, xpx_energybounds, tau);
    }
    else if (argc == 2){ 
    // Obtain trajectory and function M for initial conditions from a file
        string filename_initconds = "";
        double parameters[5], M_forw, M_back; // init_conds[4]
        string line;
        vector<double> time_interval(2);

        filename_initconds = argv[1];
        cout << "Reading file of initial conditions: " << filename_initconds << endl;

        ifstream readfile;
        readfile.open(filename_initconds);
    
        if (!readfile) {
            cerr << "Unable to open file";
            exit(GSL_FAILURE);   // call system to stop
        }
        
        // Read the parameters
        getline(readfile, line);
        istringstream param_data(line);
        // Parameters in this order: number of initial conditions, initial time, final time, forward integration flag, backward integration flag, total energy
        param_data >> parameters[0] >> parameters[1] >> parameters[2] >> parameters[3] >> parameters[4] >> parameters[5];

        // cout << parameters[1] << endl;

        vector<vector<double> > init_conds(int(parameters[0]), vector<double>(4));
        for( int i = 0; i < parameters[0] ; i++){
            getline(readfile, line);
            istringstream param_data(line);
            param_data >> init_conds[i][0] >> init_conds[i][1] >> init_conds[i][2] >> init_conds[i][3];
            
            // cout <<  std::fixed << std::setprecision(15) << setw(15) << init_conds[i][3] << endl;
        }
        readfile.close();
    

        // Call trajectory integration function 
        for (int i = 0; i < parameters[0]; i++ ){

            time_interval[0] = parameters[1];
            time_interval[1] = parameters[1] + parameters[2];

            int ic_id = i;
            if (parameters[3] == 1){
                M_forw = get_traj_funcM(time_interval, init_conds[i], ic_id);
            }
            else
            {
                M_forw = 0;
            }
            
            if (parameters[4] == 1){
                M_back = get_traj_funcM(time_interval, init_conds[i], ic_id);
            }
            else
            {
                M_back = 0;
            }
            
            double M = M_forw + M_back;

            
        }

    } 
    else if (argc == 4){
        puts("Computing Lagrangian descriptor for GIVEN resolution and integration time");
        int xRes = stoi(argv[1]);
        int yRes = stoi(argv[2]);
        double tau = stod(argv[3]);
        puts("Calculating bounds of energy surface's intersection on the (x,px) section ...");
        
        vector<double> xpx_energybounds = get_energybounds_sosatyw(E);
        // cout << xpx_energybounds[3] << endl;

        int succ  = funcM_2dof(xRes, yRes, xpx_energybounds, tau);

        // int succ  = funcM_2dof(xRes, yRes, xpx_energybounds, tau);
    }    
    else {
        puts("Incorrect number of inputs, exiting program");
        exit(GSL_FAILURE);
    }     


    return GSL_SUCCESS;
}



double quadratic (double x, void *params)
{
    struct quadratic_params *p = (struct quadratic_params *) params;

    double E = p->E;
    double D_X = p->D_X; 
    double ZETA = p->ZETA;
    double LAMBDA = p->LAMBDA;
    double EPSILON_S = p->EPSILON_S;

    return (E - (D_X*pow((1 - exp(-LAMBDA * x)), 2.0) - exp(-ZETA * LAMBDA * x) + EPSILON_S));
}

double quadratic_deriv (double x, void *params)
{
    struct quadratic_params *p = (struct quadratic_params *) params;

    double E = p->E;
    double D_X = p->D_X; 
    double ZETA = p->ZETA;
    double LAMBDA = p->LAMBDA;
    double EPSILON_S = p->EPSILON_S;

    return -(2.0*LAMBDA*D_X*exp(-LAMBDA * x)*(1 - exp(-LAMBDA * x)) + exp(-ZETA * LAMBDA * x));
}

void quadratic_fdf (double x, void *params, double *y, double *dy)
{
    struct quadratic_params *p
        = (struct quadratic_params *) params;

    double E = p->E;
    double D_X = p->D_X; 
    double ZETA = p->ZETA;
    double LAMBDA = p->LAMBDA;
    double EPSILON_S = p->EPSILON_S;

    *y = (E - (D_X*pow((1 - exp(-LAMBDA * x)), 2.0) - exp(-ZETA * LAMBDA * x) + EPSILON_S));
    *dy = -(2.0*LAMBDA*D_X*exp(-LAMBDA * x)*(1 - exp(-LAMBDA * x)) + exp(-ZETA * LAMBDA * x));
}


double Hderiv (double x, void *params)
{
    struct quadratic_params *p = (struct quadratic_params *) params;

    double E = p->E;
    double D_X = p->D_X; 
    double ZETA = p->ZETA;
    double LAMBDA = p->LAMBDA;
    double EPSILON_S = p->EPSILON_S;

    return (2.0*LAMBDA*D_X*exp(-LAMBDA * x)*(1 - exp(-LAMBDA * x)) + exp(-ZETA * LAMBDA * x));
}


double Hderiv_deriv (double x, void *params)
{
    struct quadratic_params *p = (struct quadratic_params *) params;

    double E = p->E;
    double D_X = p->D_X; 
    double ZETA = p->ZETA;
    double LAMBDA = p->LAMBDA;
    double EPSILON_S = p->EPSILON_S;

    return (2.0*pow(LAMBDA,2.0)*D_X*exp(-LAMBDA * x)*(exp(-LAMBDA * x) - LAMBDA*(1 - exp(-LAMBDA * x))) - pow((ZETA * LAMBDA),2.0)*exp(-ZETA * LAMBDA * x));
}

void Hderiv_fdf (double x, void *params, double *y, double *dy)
{
    struct quadratic_params *p
        = (struct quadratic_params *) params;

    double E = p->E;
    double D_X = p->D_X; 
    double ZETA = p->ZETA;
    double LAMBDA = p->LAMBDA;
    double EPSILON_S = p->EPSILON_S;

    *y = (2.0*LAMBDA*D_X*exp(-LAMBDA * x)*(1 - exp(-LAMBDA * x)) + exp(-ZETA * LAMBDA * x));
    *dy = (2.0*pow(LAMBDA,2.0)*D_X*exp(-LAMBDA * x)*(exp(-LAMBDA * x) - LAMBDA*(1 - exp(-LAMBDA * x))) - pow((ZETA * LAMBDA),2.0)*exp(-ZETA * LAMBDA * x));
}



/*
Energy surface's intersection with the SOS at y = yw 
*/
vector<double> get_energybounds_sosatyw(double E){

    double xMin = get_xbound_sosatyw(E, -1.0);
    double xMax = get_xbound_sosatyw(E, 0.5);
    printf("x-bounds %15.12f \t %15.12f \n", xMin, xMax);

    double pxMin = 0, pxMax = 0;

    vector<double> xpx_bounds(4);

    // Now solve for the maximum px-coordinate
    pxMax = get_pxmax_sosatyw(E, 0.0);
    pxMin = -pxMax; // due to symmetry about px = 0 axis
    printf("px-bounds %15.12f \t %15.12f \n", pxMin, pxMax);

    xpx_bounds[0] = xMin;
    xpx_bounds[1] = xMax;
    xpx_bounds[2] = pxMin;
    xpx_bounds[3] = pxMax;

    return xpx_bounds;
}


double get_pxmax_sosatyw(double E, double guess_xcoord){

    int status;
    int iter = 0, max_iter = 10;
    const gsl_root_fdfsolver_type *T;
    gsl_root_fdfsolver *s;
    double x0, x = guess_xcoord;
    gsl_function_fdf FDF;
    struct quadratic_params params = {E, D_X, ZETA, LAMBDA, EPSILON_S};

    FDF.f = &Hderiv;
    FDF.df = &Hderiv_deriv;
    FDF.fdf = &Hderiv_fdf;
    FDF.params = &params;

    T = gsl_root_fdfsolver_newton;
    s = gsl_root_fdfsolver_alloc (T);
    gsl_root_fdfsolver_set (s, &FDF, x);

    
    // printf ("using %s method\n", gsl_root_fdfsolver_name (s));

    // printf ("%-5s %10s\n","iter", "root");

    do
        {
        status = gsl_root_fdfsolver_iterate (s);
        x0 = x;
        x = gsl_root_fdfsolver_root (s);
        status = gsl_root_test_delta (x, x0, 0, 1e-3);

        if (status == GSL_SUCCESS){
            // printf ("Converged:\n");
            // printf("x coordinate that solves dHdx(x,yw, 0,0) for guess %5.2f is %15.12f \n", guess_xcoord, x);
        }

        // printf ("%5d %15.12f \n", iter, x);
        iter++;
        }
        while (status == GSL_CONTINUE && iter < max_iter);
    

    gsl_root_fdfsolver_free (s);


    double pxMax = sqrt(2*MASS_A*(E - (D_X*pow((1 - exp(-LAMBDA * x)), 2.0) - exp(-ZETA * LAMBDA * x) + EPSILON_S)));

    return pxMax;  
}

double get_xbound_sosatyw(double E, double guess_xcoord){

    int status;
    int iter = 0, max_iter = 10;
    const gsl_root_fdfsolver_type *T;
    gsl_root_fdfsolver *s;
    double x0, x = guess_xcoord;
    gsl_function_fdf FDF;
    struct quadratic_params params = {E, D_X, ZETA, LAMBDA, EPSILON_S};

    FDF.f = &quadratic;
    FDF.df = &quadratic_deriv;
    FDF.fdf = &quadratic_fdf;
    FDF.params = &params;

    T = gsl_root_fdfsolver_newton;
    s = gsl_root_fdfsolver_alloc (T);
    gsl_root_fdfsolver_set (s, &FDF, x);

    
    // printf ("using %s method\n", gsl_root_fdfsolver_name (s));

    // printf ("%-5s %10s\n","iter", "root");

    do
        {
        status = gsl_root_fdfsolver_iterate (s);
        x0 = x;
        x = gsl_root_fdfsolver_root (s);
        status = gsl_root_test_delta (x, x0, 0, 1e-3);

        if (status == GSL_SUCCESS){
            // printf ("Converged:\n");
            // printf("x coordinate that solves H(x,yw, 0,0) for guess %5.2f is %15.12f \n", guess_xcoord, x);
        }

        // printf ("%5d %15.12f \n", iter, x);
        iter++;
        }
        while (status == GSL_CONTINUE && iter < max_iter);
    

    gsl_root_fdfsolver_free (s);

    return x;  
}


/*
Obtain one of the momenta coordinate based on given values of other coordinates on the constant energy surface
*/
double deleonberne_get_py(double x, double y, double px, double E){


    double MASS_bath = 1.0;
    double ETA1 = 1.0, ETA2 = 1.0;
    int NB = 2;
    double OMEGAC = 1.0;
    double OMEGA1 = OMEGAC*log((1 - 0.5)/NB); 
    double CX1 = sqrt((2.0*ETA1*OMEGAC)/(M_PI*NB))*OMEGA1;
    double CY1 = sqrt((2.0*ETA2*OMEGAC)/(M_PI*NB))*OMEGA1; 

    double py = 0;
    double deleonberne_pot = D_X*pow((1 - exp(-LAMBDA*x)),2.0) + (Vdagger/pow(YW,4.0))*pow(y,2.0)*(pow(y,2.0) - 2.0*pow(YW,2.0))*exp(-ZETA*LAMBDA*x) + EPSILON_S;

    if (E >= (deleonberne_pot + (1/(2.0*MASS_A))*pow(px,2.0)))
        py = -sqrt( 2.0*MASS_B*(E - (deleonberne_pot + (1/(2.0*MASS_A))*pow(px,2.0)) ) );
    else{
        py = NAN;
        // cout << "vy isn't real!" j<< endl;
    }
    return py;

}
   

/*
* Lagrangian descriptors for 2 DOF system (4D phase space) 
*/
int funcM_2dof(int xRes, int yRes, vector<double> xpx_energybounds, double tau){

    // Fatten the domain to compute LD
    double xMin = xpx_energybounds[0] + 0.01*xpx_energybounds[0]; 
    double xMax = xpx_energybounds[1] + 0.01*xpx_energybounds[1];
    double yMin = xpx_energybounds[2] + 0.01*xpx_energybounds[2];
    double yMax = xpx_energybounds[3] + 0.01*xpx_energybounds[3];;
    double ySlice = YW;
    cout << xMin << " " << xMax << " " << yMin << " " << yMax << endl;

    int ic_id;

    double q1, q2, p1, p2; // 2-DOF
    double lam = 1;

    // Define domain parameters: (x,y) -> (x, px)
    // int xRes = 100;
    // int yRes = 100;
    double delta_x = (xMax - xMin)/xRes;
    double delta_y = (yMax - yMin)/yRes;
    // double delta_x = 0.0025;
    // double delta_y = 0.0025; 
    // int xRes = (xMax - xMin)/delta_x;
    // int yRes = (yMax - yMin)/delta_y;
    // cout << delta_x << " " << delta_y << endl;
    cout << "Total initial conditions: "<< (xRes + 1)*(yRes + 1) << endl;


    // Function M integration parameters
    double t0 = 0;
    // double tau = 50.0;
    double delta_t = 1e-2; // change this time step for speed and for more accurate event location and escape time 
    // int dim = 2

    vector<double> time_interval(2);
    
    cout << "Integration time span: " << tau << endl;
    // 2D array to store grid resolution and function M
    vector<vector<double> > xMesh(yRes+1, vector<double>(xRes+1));
    vector<vector<double> > yMesh(yRes+1, vector<double>(xRes+1));
    vector<vector<double> > forwM(yRes+1, vector<double>(xRes+1));
    vector<vector<double> > backM(yRes+1, vector<double>(xRes+1));
    // vector<vector<double> > M(yRes+1, vector<double>(xRes+1));
    // vector<vector<double> > Et(yRes+1, vector<double>(xRes+1));

    
    vector<double> init_cond(dim);

    // Writing domain and LD parameters
    ofstream params_dataOut;
    params_dataOut.open("params_deleonberne_M" + to_string(yRes) + "x" + to_string(xRes) + "_finalT" + to_string(float(t0 + tau)) + ".txt", ofstream::out | ofstream::trunc);

    params_dataOut << std::fixed << std::setprecision(5) << setw(5) << xMin << "\t" << xMax << "\t" << delta_x << "\t" << yMin << "\t" << yMax << "\t" << delta_y << "\t" << t0 << "\t" << tau << "\t" << delta_t << "\t" << E << "\n";

    params_dataOut.close();

    
    clock_t timer;
    timer = clock();
    for (int i = 0; i < yRes + 1; i++){
        for (int j = 0; j < xRes + 1; j++){
            
            ic_id = i*(xRes + 1) + j;
            // cout << ic_id << endl;
            // cout << init_cond[0] << "\t" << init_cond[2] << endl;

            xMesh[i][j] = init_cond[0];
            yMesh[i][j] = init_cond[2];

            // 2-DOF: Generate the other two coordinates using the constraint that the LD is being computed on a section and energy is conserved
            init_cond[0] = xMin + j*delta_x;
            init_cond[1] = ySlice;
            init_cond[2] = 0.0;
            init_cond[3] = 0.0;
            init_cond[4] = yMin + i*delta_y;
            // init_cond[5] = py;
            init_cond[6] = 0.0;
            init_cond[7] = 0.0;
            init_cond[8] = 0.0; // initialize the value of 
            double momenta = deleonberne_get_py(init_cond[0], init_cond[1], init_cond[2], E);
            if (isnan(momenta)){
                forwM[i][j] = NAN;
                backM[i][j] = NAN;
                // M[i][j] = NAN;
                // Et[i][j] = NAN;
                // puts("Not on the energy surface");
            }
            else {
                init_cond[5] = momenta;
            
                // Solve for forward trajectory 
                // Integrate from t_0 to t_0 + tau to obtain forward function M
                double M_forw = 0;
                double escape_timeF = t0 + tau;
                time_interval[0] = t0;
                time_interval[1] = t0 + tau;
                M_forw = gsl_integrate(time_interval, init_cond, delta_t, ic_id);
                // gsl_LD_escape(time_interval, init_cond, delta_t, ic_id, M_forw, escape_timeF);

                // Solve for backward trajectory
                // Integrate from t_0 to t_0 - tau to obtain backward function M
                double M_back = 0;
                double escape_timeB = t0 - tau;
                time_interval[0] = t0;
                time_interval[1] = t0 - tau;
                M_back = gsl_integrate(time_interval, init_cond, -delta_t, ic_id);
                // gsl_LD_escape(time_interval, init_cond, -delta_t, ic_id, M_back, escape_timeB);

                // cout << M_forw << " " << M_back << endl;
                // cout << escape_timeB << " " << escape_timeF << endl;

                // Save the forward and backward separately
                forwM[i][j] = M_forw;
                backM[i][j] = M_back;

                // Sum forward and backward and save as an array
                // M[i][j] = M_forw + M_back;
                // Et[i][j] = escape_timeF + fabs(escape_timeB);
            }

        }
    }

    timer = clock() - timer;
    printf ("It took me %lf secs.\n",((double)timer)/CLOCKS_PER_SEC);


    ofstream Mf_dataOut;
    Mf_dataOut.open("deleonberne_forwM" + to_string(yRes) + "x" + to_string(xRes) + "_finalT" + to_string(float(t0 + tau)) + ".txt", ofstream::out | ofstream::trunc);
    for (int i = 0; i < yRes + 1; i++){
        for (int k = 0; k < xRes; k++)
            Mf_dataOut << std::fixed << std::setprecision(15) << setw(15) << forwM[i][k] << "\t";

        Mf_dataOut << std::fixed << std::setprecision(15) << setw(15) << forwM[i][xRes] << endl;
    }
    Mf_dataOut.close();

    ofstream Mb_dataOut;
    Mb_dataOut.open("deleonberne_backM" + to_string(yRes) + "x" + to_string(xRes) + "_finalT" + to_string(float(t0 + tau)) + ".txt", ofstream::out | ofstream::trunc);
    for (int i = 0; i < yRes + 1; i++){
        for (int k = 0; k < xRes; k++)
            Mb_dataOut << std::fixed << std::setprecision(15) << setw(15) << backM[i][k] << "\t";

        Mb_dataOut << std::fixed << std::setprecision(15) << setw(15) << backM[i][xRes] << endl;
    }
    Mb_dataOut.close();


    // Writing Lagrangian Descriptor data
    // ofstream M_dataOut;
    // M_dataOut.open("deleonberne_M" + to_string(yRes) + "x" + to_string(xRes) + "_finalT" + to_string(float(t0 + tau)) + ".txt", ofstream::out | ofstream::trunc);
    // for (int i = 0; i < yRes + 1; i++){
    //     for (int k = 0; k < xRes; k++)
    //         M_dataOut << std::fixed << std::setprecision(15) << setw(15) << M[i][k] << "\t";

    //     M_dataOut << std::fixed << std::setprecision(15) << setw(15) << M[i][xRes] << endl;
    // }
    // M_dataOut.close();

    // Writing escape time data
    // ofstream Et_dataOut;
    // Et_dataOut.open("deleonberne_escape" + to_string(yRes) + "x" + to_string(xRes) + "_finalT" + to_string(float(t0 + tau)) + ".txt", ofstream::out | ofstream::trunc);
    // for (int i = 0; i < xRes + 1; i++){
    //     for (int k = 0; k < yRes; k++)
    //         Et_dataOut << std::fixed << std::setprecision(15) << setw(15) << Et[i][k] << "\t";

    //     Et_dataOut << std::fixed << std::setprecision(15) << setw(15) << Et[i][yRes] << endl;
    // }
    // Et_dataOut.close();

    return GSL_SUCCESS;
}


double get_traj_funcM(vector<double> time_interval, vector<double> init_cond, int file_id){

    double Mx = 0;
    double tdelta = 0.01;
    double t, t_next;	/* current and next independent variable */
    double tmin, tmax;	/* integration time span, lower and upper bounds */
    
    // Time span
    tmin = time_interval[0];            /* starting t value */
    tmax = time_interval[1];			/* final t value */
    
    int time_steps = int((tmax - tmin)/tdelta);

    // cout << "Start time: " << tmin << ", End time: " << tmax << endl;
    // cout << "Start time: " << time_interval[0] << ", End time: " << time_interval[1] << endl;

    double del = ZETA;
    /* load values into the my_system structure */
    gsl_odeiv2_system sys = {deleonberne_rhs, NULL, dim, &del};

    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);
    /* define the type of routine for making steps: */
    /* some other possibilities (see GSL manual):          
       = gsl_odeiv_step_rk4;
       = gsl_odeiv_step_rkck;
       = gsl_odeiv_step_rk8pd;
       = gsl_odeiv_step_rk4imp;
       = gsl_odeiv_step_bsimp;  
       = gsl_odeiv_step_gear1;
       = gsl_odeiv_step_gear2;
    */

    clock_t timer;
    

    /* 
       allocate/initialize the stepper, the control function, and the
       evolution function.
    */
    // gsl_odeiv_step *step_ptr = gsl_odeiv_step_alloc (type_ptr, dim);
    // gsl_odeiv_control *control_ptr = gsl_odeiv_control_y_new (eps_abs, eps_rel);
    // gsl_odeiv_evolve *evolve_ptr = gsl_odeiv_evolve_alloc (dim);

    // gsl_odeiv_system my_system;	/* structure with the rhs function, etc. */

    double mu = del;		 /* parameter for the diffeq */
    vector<double> y(dim);  /* current solution vector */

    /* initial values */
    for (int k = 0; k < dim; k++){
        y[k] = init_cond[k];
        // cout << init_cond[k] << "\t";
    }
    // cout << endl;


    
    // Generate the array to store the trajectory
    vector<vector<double> > trajectory(time_steps + 1, vector<double>(dim + 1));

    t = tmin; /* initialize t */
    // printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);	/* initial values */
    int k = 0; // time step counter
    
    
    trajectory[k][0] = t;
    for (int dim_cnt = 1; dim_cnt < dim + 1; dim_cnt++ ){
        trajectory[k][dim_cnt] = y[dim_cnt-1];
        // printf("%.5e \t", trajectory[k][dim_cnt]);     
    }

    // puts("Begin trajectory computation ... "); 
    timer = clock();
    
    /* step to tmax from tmin */
    t_next = tmin + tdelta;
    // while ( t_next <= tmax )
    while ( ( k < time_steps ) & (fabs(y[0]) < 1e3) )
    // for (t_next = tmin + delta_t; t_next <= tmax; t_next += delta_t)
    {

        // printf("%.5e \t %.5e \t %d \n", t, t_next, k + 1);
        // double initial_ev_val = eventfun(y); 
        double t_curr = t;
        vector<double> y_curr(dim);
        y_curr = y;
        // copy( y, y+dim, y_curr.begin() );
        
        double *y_start = &y[0]; // because y is STL vector

        while (t < t_next || t > t_next)	/* evolve from t to t_next */
        {
            int status = gsl_odeiv2_driver_apply (d, &t, t_next, y_start);


        if (status != GSL_SUCCESS)
        {
            // printf ("error, return value=%d\n", status);
            break;
            // exit;
        }

        }
        
        // printf ("%.5e \t %.5e \t%.5e\n", t, y[0], y[1]); /* print at t=t_next */
        // printf ("%.5e \t %.5e \t%.5e\n", t, y_start[0], y_start[1]); /* print at t=t_next */

        // if ( time_steps > 1000 & k % 1000 == 0 )
            // printf("Number of time steps completed: %d \n", k);

        k = k + 1;
        t_next = t_next + tdelta;
    
        trajectory[k][0] = t;
        for (int dim_cnt = 1; dim_cnt < dim + 1; dim_cnt++ ){
            trajectory[k][dim_cnt] = y_start[dim_cnt-1];
        }
        // Mx += trajectory[k][dim];

        // Saving trajectory to file
        int succ = write_matrix(trajectory, tmax, mu, file_id);

        fflush(stdout);
    } // end of trajectory generation loop

    timer = clock() - timer;
    // printf ("It took me %lf secs.\n",((double)timer)/CLOCKS_PER_SEC);

    /* all done; free up the gsl_odeiv stuff */
    gsl_odeiv2_driver_free (d);

    Mx = abs(trajectory[k][dim]);
    cout << "Value of M: " << Mx << endl;

    return Mx; 

}



/*********** INTEGRATING USING GSL ODE SOLVERS *****************/ 
double gsl_integrate(vector<double> time_interval, vector<double> init_cond, double delta_t, int file_id){

    // Event catching parameters
    // int ev_dir =  1; //-1 for detecting positive->negative 0-crossing, +1 for detecting negative->positive 0-crossing

    double Mx = 0;
    double t, t_next;	/* current and next independent variable */
    double tmin, tmax;	/* integration time span, lower and upper bounds */
    
    // Time span
    tmin = time_interval[0];            /* starting t value */
    tmax = time_interval[1];			/* final t value */
    
    // cout << "Start time: " << tmin << ", End time: " << tmax << endl;
    // cout << "Start time: " << time_interval[0] << ", End time: " << time_interval[1] << endl;

    // double delta_t = 1e-4; // time step at which the trajectory is saved
    int time_steps = (tmax - tmin)/delta_t;

    // int N = pow(2, 5); // the maximum number of segments 
    // int time_steps = floor(log2(N)) + 1;
    // double delta_t = (tmax - tmin)/time_steps;
    // printf("Time steps: %d\n", time_steps);

    double h = (delta_t/abs(delta_t))*1e-6; /* starting step size for ode solver */

    // int dim = 3; // order of ODEs + 1, derive this from size of initial condition vector 
    double eps_abs = 1.e-8;	    /* absolute error requested */
    double eps_rel = 1.e-10;	/* relative error requested */
    double lam = 1.0;

    clock_t timer;
    /* define the type of routine for making steps: */
    const gsl_odeiv_step_type *type_ptr = gsl_odeiv_step_rkf45;
    /* some other possibilities (see GSL manual):          
       = gsl_odeiv_step_rk4;
       = gsl_odeiv_step_rkck;
       = gsl_odeiv_step_rk8pd;
       = gsl_odeiv_step_rk4imp;
       = gsl_odeiv_step_bsimp;  
       = gsl_odeiv_step_gear1;
       = gsl_odeiv_step_gear2;
    */

    /* 
       allocate/initialize the stepper, the control function, and the
       evolution function.
    */
    gsl_odeiv_step *step_ptr = gsl_odeiv_step_alloc (type_ptr, dim);
    gsl_odeiv_control *control_ptr = gsl_odeiv_control_y_new (eps_abs, eps_rel);
    gsl_odeiv_evolve *evolve_ptr = gsl_odeiv_evolve_alloc (dim);

    gsl_odeiv_system my_system;	/* structure with the rhs function, etc. */

    double mu = lam;		        /* parameter for the diffeq */
    // double y[dim];			/* current solution vector */
    vector<double> y(dim);
    // cout << dim << endl;
    // cout << init_cond.size() << endl;

    /* initial values */
    for (int k = 0; k < dim; k++){
        y[k] = init_cond[k];
        // cout << init_cond[k] << "\t";
    }
    // cout << endl;

    /* load values into the my_system structure */
    // my_system.function = auto_rot_saddle_pt;	/* the right-hand-side functions dy[i]/dt */
    my_system.function = deleonberne_rhs;
    my_system.jacobian = NULL;	/* the Jacobian df[i]/dy[j] */
    my_system.dimension = dim;	/* number of diffeq's */
    my_system.params = &mu;	    /* parameters to pass to rhs and jacobian */
    
    // Generate the array to store the trajectory
    vector<vector<double> > trajectory(time_steps + 1, vector<double>(dim + 1));

    t = tmin; /* initialize t */
    // printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);	/* initial values */

    
    int k = 0; // time step counter
    trajectory[k][0] = t;
    for (int dim_cnt = 1; dim_cnt < dim + 1; dim_cnt++ ){
        trajectory[k][dim_cnt] = y[dim_cnt-1];
        // printf("%.5e \t", trajectory[k][dim_cnt]);     
    }

    // 2D array to store events
    // vector<vector<double> > event_states;

    // puts("Begin trajectory computation ... "); 
    timer = clock();
    int l = 0; // event counter
    /* step to tmax from tmin */
    t_next = tmin + delta_t;
    // while ( t_next <= tmax )
    while ( k < time_steps )
    // for (t_next = tmin + delta_t; t_next <= tmax; t_next += delta_t)
    {

        // printf("%.5e \t %.5e \t %d \n", t, t_next, k + 1);
        // double initial_ev_val = eventfun(y); 
        double t_curr = t;
        vector<double> y_curr(dim);
        y_curr = y;
        // copy( y, y+dim, y_curr.begin() );
        
        double *y_start = &y[0]; // because y is STL vector

        while (t < t_next || t > t_next)	/* evolve from t to t_next */
        {
            gsl_odeiv_evolve_apply (evolve_ptr, control_ptr, step_ptr,
                                    &my_system, &t, t_next, &h, y_start);
        }
        
        // printf ("%.5e \t %.5e \t%.5e\n", t, y[0], y[1]); /* print at t=t_next */
        // printf ("%.5e \t %.5e \t%.5e\n", t, y_start[0], y_start[1]); /* print at t=t_next */

        // if ( time_steps > 1000 & k % 1000 == 0 )
            // printf("Number of time steps completed: %d \n", k);

        k = k + 1;
        t_next = t_next + delta_t;
    
        trajectory[k][0] = t;
        for (int dim_cnt = 1; dim_cnt < dim + 1; dim_cnt++ ){
            trajectory[k][dim_cnt] = y_start[dim_cnt-1];
        }
        // Mx += trajectory[k][dim];

        // Saving trajectory to file
        // int succ = write_matrix(trajectory, tmax, mu, file_id);

        fflush(stdout);
    } // end of trajectory generation loop

    timer = clock() - timer;
    // printf ("It took me %lf secs.\n",((double)timer)/CLOCKS_PER_SEC);

    /* all done; free up the gsl_odeiv stuff */
    gsl_odeiv_evolve_free (evolve_ptr);
    gsl_odeiv_control_free (control_ptr);
    gsl_odeiv_step_free (step_ptr);

    // Mx = M_Lopesino2017(trajectory, lam, time_steps);

    Mx = abs(trajectory[k][dim]);
    // cout << "Value of M: " << Mx << endl;

    return Mx;    
}


/* 
Integration using GSL solvers along with function M and stopping when escape occurs by leaving a prescribed domain in configuration space
*/

int gsl_LD_escape(vector<double> time_interval, vector<double> init_cond, double delta_t, int file_id, double& Mt, double& escape_time){

    // Event catching parameters
    // int ev_dir =  1; //-1 for detecting positive->negative 0-crossing, +1 for detecting negative->positive 0-crossing

    // double Mx = 0;
    double t, t_next;	/* current and next independent variable */
    double tmin, tmax;	/* integration time span, lower and upper bounds */
    
    // Time span
    tmin = time_interval[0];            /* starting t value */
    tmax = time_interval[1];			/* final t value */
    
    // cout << "Start time: " << tmin << ", End time: " << tmax << endl;
    // cout << "Start time: " << time_interval[0] << ", End time: " << time_interval[1] << endl;

    // double delta_t = 1e-4; // time step at which the trajectory is saved
    int time_steps = (tmax - tmin)/delta_t;

    // int N = pow(2, 5); // the maximum number of segments 
    // int time_steps = floor(log2(N)) + 1;
    // double delta_t = (tmax - tmin)/time_steps;
    // printf("Time steps: %d\n", time_steps);

    double h = (delta_t/abs(delta_t))*1e-6; /* starting step size for ode solver */

    // int dim = 3; // order of ODEs + 1, derive this from size of initial condition vector 
    double eps_abs = 1.e-8;	    /* absolute error requested */
    double eps_rel = 1.e-10;	/* relative error requested */
    double lam = 1.0;

    clock_t timer;
    /* define the type of routine for making steps: */
    const gsl_odeiv_step_type *type_ptr = gsl_odeiv_step_rkf45;
    /* some other possibilities (see GSL manual):          
       = gsl_odeiv_step_rk4;
       = gsl_odeiv_step_rkck;
       = gsl_odeiv_step_rk8pd;
       = gsl_odeiv_step_rk4imp;
       = gsl_odeiv_step_bsimp;  
       = gsl_odeiv_step_gear1;
       = gsl_odeiv_step_gear2;
    */

    /* 
       allocate/initialize the stepper, the control function, and the
       evolution function.
    */
    gsl_odeiv_step *step_ptr = gsl_odeiv_step_alloc (type_ptr, dim);
    gsl_odeiv_control *control_ptr = gsl_odeiv_control_y_new (eps_abs, eps_rel);
    gsl_odeiv_evolve *evolve_ptr = gsl_odeiv_evolve_alloc (dim);

    gsl_odeiv_system my_system;	/* structure with the rhs function, etc. */

    double mu = lam;		        /* parameter for the diffeq */
    // double y[dim];			/* current solution vector */
    vector<double> y(dim);
    // cout << dim << endl;
    // cout << init_cond.size() << endl;

    /* initial values */
    for (int k = 0; k < dim; k++){
        y[k] = init_cond[k];
        // cout << init_cond[k] << "\t";
    }
    // cout << endl;

    /* load values into the my_system structure */
    // my_system.function = auto_rot_saddle_pt;	/* the right-hand-side functions dy[i]/dt */
    my_system.function = deleonberne_rhs;
    my_system.jacobian = NULL;	/* the Jacobian df[i]/dy[j] */
    my_system.dimension = dim;	/* number of diffeq's */
    my_system.params = &mu;	    /* parameters to pass to rhs and jacobian */
    
    // Generate the array to store the trajectory
    vector<vector<double> > trajectory(time_steps + 1, vector<double>(dim + 1));

    t = tmin; /* initialize t */
    // printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);	/* initial values */

    
    int k = 0; // time step counter
    trajectory[k][0] = t;
    for (int dim_cnt = 1; dim_cnt < dim + 1; dim_cnt++ ){
        trajectory[k][dim_cnt] = y[dim_cnt-1];
        // printf("%.5e \t", trajectory[k][dim_cnt]);     
    }

    // 2D array to store events
    // vector<vector<double> > event_states;

    // puts("Begin trajectory computation ... "); 
    timer = clock();
    int l = 0; // event counter
    /* step to tmax from tmin */
    t_next = tmin + delta_t;
    // while ( t_next <= tmax )
    while ( k < time_steps )
    // for (t_next = tmin + delta_t; t_next <= tmax; t_next += delta_t)
    {

        // printf("%.5e \t %.5e \t %d \n", t, t_next, k + 1);
        // double initial_ev_val = eventfun(y); 
        double t_curr = t;
        vector<double> y_curr(dim);
        y_curr = y;
        // copy( y, y+dim, y_curr.begin() );
        
        double *y_start = &y[0]; // because y is STL vector

        while (t < t_next || t > t_next)	/* evolve from t to t_next */
        {
            gsl_odeiv_evolve_apply (evolve_ptr, control_ptr, step_ptr,
                                    &my_system, &t, t_next, &h, y_start);
        }
        
        // printf ("%.5e \t %.5e \t%.5e\n", t, y[0], y[1]); /* print at t=t_next */
        // printf ("%.5e \t %.5e \t%.5e\n", t, y_start[0], y_start[1]); /* print at t=t_next */

        // if ( time_steps > 1000 & k % 1000 == 0 )
            // printf("Number of time steps completed: %d \n", k);

    
        trajectory[k][0] = t;
        for (int dim_cnt = 1; dim_cnt < dim + 1; dim_cnt++ ){
            trajectory[k][dim_cnt] = y_start[dim_cnt-1];
        }
        // Mt += trajectory[k][dim];

        Mt = abs(trajectory[k][dim]);
        escape_time = t_next;

        // Saving trajectory to file
        // int succ = write_matrix(trajectory, tmax, mu, file_id);

        if (y_start[1] < -8.0 || y_start[1] > -6.0 || y_start[0] > 6.5 || y_start[0] < 4.5) {
            // puts("Escape event has occurred.");
            break;
        } // E = 15.250

        
        // if (y_start[1] < -10.0 || y_start[1] > -5.0 || y_start[0] > 8.0 || y_start[0] < 3.0) {
        //     // puts("Escape event has occurred.");
        //     break;
        // } // E = 16.250

        // if (fabs(y_start[1]) > 10e3 || fabs(y_start[0]) > 10e3) {
        //     // puts("Escape event has occurred.");
        //     break;
        // } // Escape to infinity
        
        k = k + 1;
        t_next = t_next + delta_t;
        
        fflush(stdout);
    } // end of trajectory generation loop

    timer = clock() - timer;
    // printf ("It took me %lf secs.\n",((double)timer)/CLOCKS_PER_SEC);

    /* all done; free up the gsl_odeiv stuff */
    gsl_odeiv_evolve_free (evolve_ptr);
    gsl_odeiv_control_free (control_ptr);
    gsl_odeiv_step_free (step_ptr);

    
    // cout << "Value of M: " << Mt << endl;
    // cout << "Escape time: " << fabs(escape_time) << endl;

    return GSL_SUCCESS;    
}

/*==================== Writing matrix to file ==============================*/
int write_matrix(const vector<vector<double> >& trajectory, double tmax, double mu, int file_id){

    int time_steps = trajectory.size();
    int num_cols = trajectory[0].size(); 
    // cout << time_steps << " " << num_cols << endl;

    // Writing trajectory data
    ofstream trajDataOut;
    trajDataOut.open("test_traj" + to_string(int(file_id)) + "_finalT" + to_string(int(tmax)) + "_lam"+ to_string(int(mu)) + ".txt", ofstream::out | ofstream::trunc);
    for (int i = 0; i < time_steps; i++){
        for (int k = 0; k < num_cols; k++)
            // trajDataOut << std::fixed << std::setprecision(15) << setw(15) << trajectory[i][k] << "\t";
            trajDataOut << std::scientific << std::setprecision(15) << setw(15) << trajectory[i][k] << "\t";
        trajDataOut << endl;
    }
    trajDataOut.close();

    return GSL_SUCCESS;

}


/*
    * Velocity field for dynamics governed by Barbanis potential 
    * Phase space is 4D, with x, y, px, py
*/
int deleonberne_rhs(double t, const double posIn[], double velOut[], void *params){

    // double De = 100, k = 200, re = 1;
    double b = *(double *) params;
    double x, y, px, py;
    double dV_dx, dV_dy;
    double p = 0.5;

    x = posIn[0];
    y = posIn[1];

    dV_dx = -( 
        2*D_X*LAMBDA*exp(-LAMBDA*x)*(exp(-LAMBDA*x) - 1) + 
        (Vdagger/pow(YW,4.0))*ZETA*LAMBDA*pow(y,2.0)*(pow(y,2.0) - 2.0*pow(YW,2.0))*exp(-ZETA*LAMBDA*x) );

    dV_dy = ( 4*((Vdagger/pow(YW,4.0))*y*( pow(y,2.0) - pow(YW,2.0) )*exp(-ZETA*LAMBDA*x)) );

    velOut[0] = (posIn[2]/MASS_A);
    velOut[1] = (posIn[3]/MASS_B);
    velOut[2] = -(dV_dx);
    velOut[3] = -(dV_dy);
    velOut[4] = 0;
    for (int i = 0; i < dim - 1; i++){
        velOut[4] += pow(abs(velOut[i]),p);
    }

    // printf("%f \t %f \t %f \t %f\n", velOut[0], velOut[1], velOut[2], velOut[3]);

    return GSL_SUCCESS;

}


/*
    * Velocity field for dynamics governed by De Leon-Berne system in a bath
    * Phase space is 8 dimensional
*/
int deleonberne_bath_rhs(double t, const double posIn[], double velOut[], void *params){

    // double De = 100, k = 200, re = 1;
    double MASS_bath = 1.0;
    double ETA1 = 1.0, ETA2 = 1.0;
    int NB = 2;
    double OMEGAC = 1.0;
    double OMEGA1 = OMEGAC*log((1 - 0.5)/NB); 
    double CX1 = sqrt((2.0*ETA1*OMEGAC)/(M_PI*NB))*OMEGA1;
    double CY1 = sqrt((2.0*ETA2*OMEGAC)/(M_PI*NB))*OMEGA1; 

    double b = *(double *) params;
    double x, y, x1, y1, px, py, px1, py1;
    double dV_dx, dV_dy, dV_dxj, dV_dyj;
    double p = 0.5;

    x = posIn[0];
    y = posIn[1];
    x1 = posIn[2];
    y1 = posIn[3];

    dV_dx = -( 
        2*D_X*LAMBDA*exp(-LAMBDA*x)*(exp(-LAMBDA*x) - 1) + 
        (Vdagger/pow(YW,4.0))*ZETA*LAMBDA*pow(y,2.0)*(pow(y,2.0) - 2.0*pow(YW,2.0))*exp(-ZETA*LAMBDA*x) + (CX1/OMEGA1)*( OMEGA1*x1 - (CX1*x)/OMEGA1 ) 
        );

    dV_dy = ( 
        4*((Vdagger/pow(YW,4.0))*y*( pow(y,2.0) - pow(YW,2.0) )*exp(-ZETA*LAMBDA*x)) - 
        (CY1/OMEGA1)*( OMEGA1*y1 - (CY1*y)/OMEGA1 )
        );

    dV_dxj = OMEGA1*( OMEGA1*x1 - (CX1*x)/OMEGA1 );

    dV_dyj = OMEGA1*( OMEGA1*y1 - (CY1*y)/OMEGA1 );


    velOut[0] = (posIn[4]/MASS_A);
    velOut[1] = (posIn[5]/MASS_B);
    velOut[2] = (posIn[5]/MASS_bath); 
    velOut[3] = (posIn[5]/MASS_bath);
    velOut[4] = -(dV_dx);
    velOut[5] = -(dV_dy);
    velOut[6] = -(dV_dxj);
    velOut[7] = -(dV_dyj); 
    velOut[8] = 0;
    for (int i = 0; i < dim - 1; i++){
        velOut[8] += pow(abs(velOut[i]),p);
    }

    // printf("%f \t %f \t %f \t %f\n", velOut[0], velOut[1], velOut[2], velOut[3]);

    return GSL_SUCCESS;

}












