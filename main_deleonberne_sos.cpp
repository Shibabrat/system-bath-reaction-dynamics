/* UNDER DEVELOPMENT
@Purpose    C++ code to compute surface-of-section of trajectories for 
            De Leon-Berne model of isomerization

@author Shibabrat Naik [2019]
*/

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_odeiv2.h>

#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm> 
// #include <sstream>

using namespace std;

size_t dim = 4;
// Parameters
double OMEGA_X = 1.0, OMEGA_Y = 1.1, DELTA = -0.11, M = 1.0; 
double E = 11.000;
double lambda = 1.0, omega2 = 1.0;


/* Function prototypes */
double get_py_linear_symp_2dof(double q1, double p2, double E);
double barbanis_get_vy(double x, double y, double vx, double E);
int gsl_integrate(vector<double> parameters);
int locate_event(const double ti, const double tf, vector<double> yi, vector<double> yf, double& te, vector<double>& ye);
double eventfun(vector<double> Y);
int barbanis_rhs(double t, const double posIn[], double velOut[], void *params);
int linear_symp_2dof (double t, const double posIn[], double velOut[], void *params_ptr);


// Functions not being used
int barbanis_sos();
void vel_verlet_evolve(const double h, const vector<double> yi, vector<double>& yf);
void locate_event_vv(const double ti, const double tf, vector<double> yi, vector<double> yf, double& te, vector<double>& ye);
// Functions not being used


int main(int argc, char *argv[]){

    if (argc == 1){
        puts("Incorrect number of inputs, exiting program");
        exit(GSL_FAILURE);
    }
    
    if (argc == 2){

        // Obtain trajectory and function M for initial conditions from a file
        string filename_TE_boundary = "";
        double M_forw, M_back; // init_conds[4]
        string line;
        vector<double> parameters(5);
        
        filename_TE_boundary = argv[1];
        cout << "Reading file of energy boundary: " << filename_TE_boundary << endl;

        ifstream readfile;
        readfile.open(filename_TE_boundary);

        if (!readfile) {
            cerr << "Unable to open file";
            exit(GSL_FAILURE);   // call system to stop
        }
        
        // Read the parameters
        getline(readfile, line);
        istringstream param_data(line);
        // Parameters in this order: total energy, minimum and maximum values of x, minimum and maximum values of px
        param_data >> parameters[0] >> parameters[1] >> parameters[2] >> parameters[3] >> parameters[4];

        cout << parameters[2] << endl;
        readfile.close();


        int succ  = gsl_integrate(parameters);
        // int succ  = barbanis_sos();

    }
    
    return GSL_SUCCESS;
}


/*
Obtain the positive y-momentum coordinate for given x, px and total energy values and at y = 0
*/
double get_py_linear_symp_2dof(double q1, double p1, double E){

    double py;

    double a = (0.5*lambda + omega2);
    double b = (-lambda*q1 + lambda*p1 + omega2*p1);
    double c = (0.5*lambda*(pow(q1,2.0) - 2*q1*p1) + 0.5*omega2*pow(p1,2.0) - E );

    double py1 = (-b + sqrt( pow(b,2.0) - 4.0*a*c ))/(2.0*a);
    double py2 = (-b - sqrt( pow(b,2.0) - 4.0*a*c ))/(2.0*a);

    if (py1 > 0)
        py = py1;
    else if (py2 > 0) 
        py = py2;
    else{
        py = NAN;
        cout << "py isn't real!" << endl;
    }

    return py;
}

/*
Obtain one of the coordinate based on constant energy
*/

double barbanis_get_vy(double x, double y, double vx, double E){

    double vy;
    double barbanis_pot = 0.5*pow(OMEGA_X,2.0)*pow(x,2.0) + 0.5*pow(OMEGA_Y,2.0)*pow(y,2.0) + DELTA*x*pow(y,2.0);

    if ( 2.0*E >=  (2.0*barbanis_pot + pow(vx,2.0)) ){
        vy = sqrt( 2.0*E - 2.0*barbanis_pot - pow(vx,2.0) );
    }
    else{
        vy = NAN;
    }
    

    return vy;
}


/* 
 * INTEGRATING USING GSL ODE SOLVERS 
 * */ 
int gsl_integrate(vector<double> parameters){

    // Event catching parameters
    int ev_dir =  1; //-1 for detecting positive->negative 0-crossing, +1 for detecting negative->positive 0-crossing

    // int dim = 4;                /* order of ODEs */
    double eps_abs = 1.e-8;	    /* absolute error requested */
    double eps_rel = 1.e-10;	/* relative error requested */

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

    double mu = DELTA;		/* parameter for the diffeq */
    double total_energy = parameters[0];   // Total energy
    // double y[dim];			/* current solution vector */
    
    // Initial condition along each coordinate
    int num_ics = 10;

    // Limits of the (x,px) for linear symplectic tranformed system
    // double xMin = -1.0;
    // double xMax = 1.0;
    // double pxMin = -1.0;
    // double pxMax = 1.0;

    // Limits of the (x,px) for E = 0.5
    // double xMin = -0.999;
    // double xMax = 0.999;
    // double pxMin = -0.999;
    // double pxMax = 0.999;

    // Limits of the (x, px) surface for E = 10.000
    double xMin = parameters[1];
    double xMax = parameters[2];
    double pxMin = parameters[3];
    double pxMax = parameters[4];

    // Limits of the (x,px) for E = 16.375
    // double xMin = -5.722;
    // double xMax = 5.722;
    // double pxMin = -5.722;
    // double pxMax = 5.722;

    // Limits of the (x,px) for E = 15.250
    // double xMin = -5.522;
    // double xMax = 5.522;
    // double pxMin = -5.522;
    // double pxMax = 5.522;

    // Limits of the (x,px) for E = 15.125
    // double xMin = -5.495;
    // double xMax = 5.495;
    // double pxMin = -5.495;
    // double pxMax = 5.495;
    
    // Limits of the (x,px) for E = 15.000
    // double xMin = -5.477;
    // double xMax = 5.477;
    // double pxMin = -5.477;
    // double pxMax = 5.477;

    // Limits of the (x,px) for E = 13.875
    // double xMin = -5.267;
    // double xMax = 5.267;
    // double pxMin = -5.267;
    // double pxMax = 5.267;

    // Limits of the (x,px) for E = 11.000
    // double xMin = -4.690;
    // double xMax = 4.690;
    // double pxMin = -4.690;
    // double pxMax = 4.690;

    double t, t_next;		/* current and next independent variable */
    double tmin, tmax, delta_t;	/* range of t and step size for output */

    double h = 1e-6;		/* starting step size for ode solver */

    /* load values into the my_system structure */
    // my_system.function = linear_symp_2dof;	/* the right-hand-side functions dy[i]/dt */
    my_system.function = barbanis_rhs;	/* the right-hand-side functions dy[i]/dt */
    my_system.jacobian = NULL;	/* the Jacobian df[i]/dy[j] */
    my_system.dimension = dim;	/* number of diffeq's */
    my_system.params = &mu;	/* parameters to pass to rhs and jacobian */

    tmin = 0.;			        /* starting t value */
    tmax = 750.0;			/* final t value */
    delta_t = 1e-3;

    // Time steps
    int time_steps = (tmax - tmin)/delta_t;
    printf("Time steps : %d \n", time_steps);

    /* Generate initial conditions on the energy surface */
    // for ( int ic_cnt_x = 0; ic_cnt_x < num_ics; ic_cnt_x++ ){
    int ic_cnt_px = 0;
    while ( ic_cnt_px < num_ics ){
        int ic_cnt_x = 0;
        while ( ic_cnt_x < num_ics ){

            vector<double> y(dim);

            // Generates initial conditions on (x,px) section and avoids the stable equilibrium at (0,0)
            y[0] = xMin + ic_cnt_x*((xMax - xMin)/num_ics) ;			
            y[1] = 0.0;			
            y[2] = pxMin + ic_cnt_px*((pxMax - pxMin)/num_ics) ;
            y[3] = barbanis_get_vy(y[0], y[1], y[2], total_energy);

            // Generates a vertical line of initial condition on (x,px) section
            // y[0] = 0.0; 			
            // y[1] = 0.0;			
            // y[2] = pxMin + ic_cnt_x*((pxMax - pxMin)/num_ics) ;;
            // y[3] = get_py_linear_symp_2dof(y[0], y[2], total_energy);

            if (y[3] == NAN)
                break;
            
            t = tmin; /* initialize t */
            cout << "Initial condition " << int(ic_cnt_x + (num_ics)*ic_cnt_px) << endl;
            printf ("%.5e %.5e %.5e\n", y[0], y[1], y[2]);	/* initial values */

            int num_steps = 0;  // Number of the time steps
            int l = 0;          // event counter

            // Generate the array to store the trajectory
            // vector<vector<double> > trajectory(time_steps + 1, vector<double>(dim+1));
            // trajectory[k][0] = t;
            // for (int dim_cnt = 1; dim_cnt < dim + 1; dim_cnt++ ){
            //     trajectory[k][dim_cnt] = y[dim_cnt-1];
            //     // printf("%.5e \t", trajectory[k][dim_cnt]);     
            // }

            // 2D array to store events
            vector<vector<double> > event_states;

            puts("Begin trajectory computation ... "); 
            timer = clock();
            

            /* step to tmax from tmin */
            t_next = tmin + delta_t;
            while ( t_next <= tmax )
            // for (t_next = tmin + delta_t; t_next <= tmax; t_next += delta_t)
            {

                // printf("%.5e \t %.5e \t %d \n", t, t_next, num_events + 1);
                double initial_ev_val = eventfun(y); 
                double t_curr = t;
                vector<double> y_curr(dim);
                y_curr = y;
                // copy( y, y+dim, y_curr.begin() );
                
                double *y_start = &y[0]; // because y is STL vector

                while (t < t_next)	/* evolve from t to t_next */
                {
                    gsl_odeiv_evolve_apply (evolve_ptr, control_ptr, step_ptr,
                                            &my_system, &t, t_next, &h, y_start);
                }
                // printf ("%.5e \t %.5e \t%.5e\n", t, y[0], y[1]); /* print at t=t_next */
                // printf ("%.5e \t %.5e \t%.5e\n", t, y_start[0], y_start[1]); /* print at t=t_next */

                if (abs(y[1]) > 10e3 || ((double)(clock() - timer)/CLOCKS_PER_SEC) > 10 ){
                    puts("Escape event has occurred or integration is taking too long!");
                    break;
                }
                double final_ev_val = eventfun(y);

                if  ( (initial_ev_val/fabs(initial_ev_val)) == (final_ev_val/fabs(final_ev_val)) ) {
                    num_steps += 1;
                    t_next += delta_t;

                    // trajectory[num_steps][0] = t;
                    // for (int dim_cnt = 1; dim_cnt < dim + 1; dim_cnt++ ){
                    //     trajectory[num_steps][dim_cnt] = y[dim_cnt-1];
                    //     // printf("%.5e \t", trajectory[num_steps][dim_cnt]);     
                    // }
                    // printf("\n");
                }
                else if ( ev_dir*(initial_ev_val/fabs(initial_ev_val)) < ev_dir*(final_ev_val/fabs(final_ev_val)) ) { // negative to positive crossing

                    // printf("Event value went from %.5e to %.5e \n", initial_ev_val, final_ev_val);
                    // printf("Event occurred between %.5e and %.5e \n", t_curr, t_next);
                    double te;
                    vector<double> ye(dim);
                    
                    int temp = locate_event(t_curr, t_next, y_curr, y, te, ye);
                    // printf("%.10f \t %.10f \t %.10f \t %.10f\n", ye[0], ye[1], ye[2], ye[3]);

                    event_states.push_back(vector<double>());
                    event_states[l].push_back(te);
                    for ( int dimcnt = 0; dimcnt < dim; dimcnt++ ){
                        event_states[l].push_back(ye[dimcnt]);
                    }
                    l += 1;

                }

                // if ( k % 10000 == 0 )
                //     printf("Number of time steps completed: %d\n", k);

                fflush(stdout);

                t_next += delta_t;
            }

            timer = clock() - timer;
            printf ("It took me %lf secs.\n",((double)timer)/CLOCKS_PER_SEC);

            // Writing trajectory data
            // ofstream trajDataOut;
            // trajDataOut.open("Barbanis_traj_finalT"+to_string(int(tmax))+"_E"+to_string(float(E))+".txt", ofstream::out | ofstream::trunc);
            // // trajDataOut.open("dm_traj_finalT"+to_string(int(tmax))+"_b"+to_string(int(mu))+".txt", ofstream::out );
            // for (int i=0; i< time_steps + 1; i = i + 100){
            //     trajDataOut << std::fixed << std::setprecision(10) << setw(10) << trajectory[i][0] << "\t";
            //     for (int k=1; k < dim + 1; k++)
            //         trajDataOut << std::fixed << std::setprecision(15) << setw(15) << trajectory[i][k] << "\t";
            //     trajDataOut << endl;
            // }
            // trajDataOut.close();


            // Writing events to file
            // printf("Size %d x %d\n ", l, event_states[0].size());
            printf("Number of events : %d\n", l);
            ofstream eventDataOut;
            eventDataOut.open("Barbanis_2dof_events_finalT"+to_string(int(tmax))+"_E"+to_string(float(total_energy))+"_id"+to_string(int(ic_cnt_x + (num_ics)*ic_cnt_px))+".txt", ofstream::out | ofstream::trunc);
            // eventDataOut.open("linear_symp_2dof_events_finalT"+to_string(int(tmax))+"_E"+to_string(float(E))+"_id"+to_string(int(ic_cnt_x))+".txt", ofstream::out | ofstream::trunc);
            for (int i=0; i< l; i++){
                for (int dim_cnt = 0; dim_cnt < dim + 1; dim_cnt++)
                    eventDataOut << std::fixed << std::setprecision(15) << setw(15) << event_states[i][dim_cnt] << "\t";
                eventDataOut << endl;
            }
            eventDataOut.close();

            ic_cnt_x += 1;
        }
    
        ic_cnt_px += 1;
    }
    
    /* all done; free up the gsl_odeiv stuff */
    gsl_odeiv_evolve_free (evolve_ptr);
    gsl_odeiv_control_free (control_ptr);
    gsl_odeiv_step_free (step_ptr);
    
    return 0;    
}

/************************** EVENT LOCATION ****************************/ 
int locate_event(const double ti, const double tf, vector<double> yi, vector<double> yf, double& te, vector<double>& ye){

    // int dim = 4;
    double eventTime;
    double h = 1e-8;

    // Bisection method parameters
    int iter = 0;
    int iterMax = 100;

    // printf("Event is between %.5f \t %.5f \n", ti, tf);
    // for ( int dim_cnt = 0; dim_cnt < dim; dim_cnt++ ){
    //     printf("%.5f \t %.5f \n", yi[dim_cnt], yf[dim_cnt]);    
    // }
    
    double prev_val = eventfun(yi);
    double curr_val = eventfun(yf);

    // cout << prev_val << " " << curr_val << endl;

    double tLeft = ti;
    double tRight = tf;
    double ev_t_guess, tInitial;
    vector<double> yMid, prevy;
    yMid = yi;

    double val_mid;
    // double E = 16.25;   // Total energy    
    double params = E;

    double eps_abs = 1.e-10;	    /* absolute error requested */
    double eps_rel = 1.e-12;	/* relative error requested */
    const gsl_odeiv2_step_type *type_ptr = gsl_odeiv2_step_rkf45;
    gsl_odeiv2_step *step_ptr = gsl_odeiv2_step_alloc (type_ptr, dim);
    gsl_odeiv2_control *control_ptr = gsl_odeiv2_control_y_new (eps_abs, eps_rel);
    gsl_odeiv2_evolve *evolve_ptr = gsl_odeiv2_evolve_alloc (dim);

    gsl_odeiv2_system sys = { linear_symp_2dof , NULL, dim, &params };
    // gsl_odeiv2_system sys = { barbanis_rhs, NULL, dim, &params };
    // gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4, h, 1e-12, 1e-12);

    double *y_start = &yi[0]; 
    while ( ( iter < iterMax )  ||  ( fabs(prev_val - curr_val) > 1e-10 ) ){

        ev_t_guess = 0.5*( tLeft + tRight );
        // printf("New guess for event time : %.10f\n", ev_t_guess);
        tInitial = tLeft;
        // double h = tRight - tLeft;

        // double *y_start = &yMid[0]; // because yMid is STL vector
        

        // printf("%f \t %f \n", tInitial, ev_t_guess);
        // printf("States at %f ... \n ", tInitial);
        // for ( int dimcnt=0; dimcnt<dim;dimcnt++){
        //     printf("%f \t", y_start[dimcnt]);
        // }
        // printf("\n");

        while (tInitial < ev_t_guess){
            gsl_odeiv2_evolve_apply (evolve_ptr, control_ptr, step_ptr,
                                    &sys, &tInitial, ev_t_guess, &h, y_start);
        }

        // int status = gsl_odeiv2_driver_apply_fixed_step (d, &tInitial, h, 1, y_start);

        // if (status != GSL_SUCCESS){
        //     printf("MAJOR ERROR \n"); 
        //     // break;
        // } 

        // printf("States at %f ... \n ", tInitial);
        // for ( int dimcnt=0; dimcnt<dim;dimcnt++){
        //     printf("%f \t", y_start[dimcnt]);
        // }
        // printf("\n");

        // event value at the mid point
        val_mid = eventfun(yi); 
        
        // printf("Event value is %f at guess %f \n", val_mid, ev_t_guess);
    
        if ( val_mid/fabs(val_mid) == curr_val/fabs(curr_val) ){ // move the curr point to mid point

            tRight = ev_t_guess;
            curr_val = val_mid;
            yf = yi;
            yi = yMid;
            double *y_start = &yi[0]; 
            // printf("Event is between %f and %f\n", tLeft, tRight);

        } 
        else if ( val_mid/fabs(val_mid) == prev_val/fabs(prev_val) ){ // move the prev point to mid point
            
            tLeft = ev_t_guess;
            prev_val = val_mid;
            yMid = yi;
            double *y_start = &yi[0];
            // printf("Event is between %f and %f .\n", tLeft, tRight);
            
        }

        // *eventTime = ev_t_guess;
       eventTime = ev_t_guess;

        iter += 1;
    }

    gsl_odeiv2_evolve_free (evolve_ptr);
    gsl_odeiv2_control_free (control_ptr);
    gsl_odeiv2_step_free (step_ptr);

    // gsl_odeiv2_driver_free (d);

    te = eventTime;
    ye = yMid;
    // printf("%f \t %f \t %f \t %f\n", ye[0], ye[1], ye[2], ye[3]);
    return GSL_SUCCESS;
}


/*************************** EVENT FUNCTION ****************************/
double eventfun(vector<double> Y) 

/**
*   @brief      This is the event function (0-crossings to be detected)
**/
{ 
    
    double result;
    result = (Y[1] - 0);

    return result; 
}


/*************************** RHS ****************************/
/* 
   Define the array of right-hand-side functions y[i] to be integrated.
  
   * params is a void pointer that is used in many GSL routines
   to pass parameters to a function
*/
/*
 * Velocity field for dynamics governed by Barbanis potential 
 * Phase space is 4D, with x, y, px, py
*/
int barbanis_rhs(double t, const double posIn[], double velOut[], void *params){

    // double De = 100, k = 200, re = 1;
    double b = *(double *) params;
    double x, y, px, py;
    double dV_dx, dV_dy;

    dV_dx = (pow(OMEGA_X,2.0)*posIn[0] + DELTA*pow(posIn[1],2.0));

    dV_dy = (pow(OMEGA_Y,2.0)*posIn[1] + 2.0*DELTA*posIn[0]*posIn[1]);

    velOut[0] = (posIn[2]/M);
    velOut[1] = (posIn[3]/M);
    velOut[2] = -(dV_dx);
    velOut[3] = -(dV_dy);

    // printf("%f \t %f \t %f \t %f\n", velOut[0], velOut[1], velOut[2], velOut[3]);

    return GSL_SUCCESS;

}


/*===== Linear symplectic transformed 2 DOF system ========================*/
int linear_symp_2dof (double t, const double posIn[], double velOut[], void *params_ptr){

    /* get parameter(s) from params_ptr; here, just a double */
    double lam = *(double *) params_ptr;

    /* evaluate the right-hand-side functions at t */
    velOut[0] = lambda*(-posIn[2]);
    velOut[1] = omega2*(+posIn[3]);
    velOut[2] = -lambda*posIn[0] + (lambda + omega2)*posIn[3];
    velOut[3] = -omega2*posIn[1] + (omega2 - lambda)*posIn[2] + 2.0*omega2*posIn[3];

    
    return GSL_SUCCESS;		/* GSL_SUCCESS defined in gsl/errno.h as 0 */
}

/*
Velocity Verlet integrator with event location for SOS 
*/

/********************* VELOCITY VERLET INTEGRATOR ************************/
// int barbanis_sos(){

//     int ev_dir =  1; //-1 for detecting positive->negative 0-crossing, +1 for detecting negative->positive 0-crossing
//     int l = 0; // event counter

//     // Time parameters
//     // unsigned long long int time_steps = 1e10;
//     double tBegin = 0;
//     double h = 1e-4;
//     double tEnd = 100;
//     unsigned long int time_steps = (int)(tEnd - tBegin)/h;
//     printf("Time steps : %d \n", time_steps);

//     clock_t timer;

//     double total_energy = 16.25; // Total energy
//     // 2D array to store events
//     vector<vector<double> > event_states;

//     vector<double> y(dim);
//     /* initial values */
//     y[0] = 1.0;			
//     y[1] = 0.0;			
//     y[2] = 0;
//     y[3] = barbanis_get_vy(y[0], y[1], y[2], total_energy);

//     double t = tBegin;
//     long int k = 0; // time step counter

//     // 2D array to store the trajectory
//     vector<vector<double> > trajectory(time_steps + 1, vector<double>(dim+1));
//     trajectory[k][0] = t;
//     for (int dim_cnt = 1; dim_cnt < dim + 1; dim_cnt++ ){
//         trajectory[k][dim_cnt] = y[dim_cnt-1];
//     }


//     printf("%.5e : \t", t);
//     for (int dim_cnt = 0; dim_cnt < dim; dim_cnt++ ){
//         printf("%.5e \t", y[dim_cnt]);     
//     }
//     printf("\n");

//     puts("Begin trajectory computation ... "); 
//     timer = clock();
//     // while ( t < tEnd ){
//     while ( k < time_steps ){

//         double initial_ev_val = eventfun(y); 
//         vector<double> y_curr = y; 
            
//         // single step using Velocity Verlet integrator
//         vel_verlet_evolve(h, y_curr, y);
        
//         // t += h;
//         // k += 1;
//         double final_ev_val = eventfun(y);

//         if  ( (initial_ev_val/fabs(initial_ev_val)) == (final_ev_val/fabs(final_ev_val)) ) {
//             t += h;
//             k += 1;
//         }
//         else if ( ev_dir*(initial_ev_val/fabs(initial_ev_val)) < ev_dir*(final_ev_val/fabs(final_ev_val)) ) { // negative to positive crossing
//             // printf("Event value went from %.5e to %.5e \n", initial_ev_val, final_ev_val);
//             // printf("Event occurred between %.5e and %.5e \n", t, t + h);
//             double te;
//             vector<double> ye(dim);
//             double t_curr = t;
//             double t_next = t + h;
            
//             // int temp = locate_event(t_curr, t_next, y_curr, y, te, ye);
//             locate_event_vv(t_curr, t_next, y_curr, y, te, ye);
//             printf("%.15f \t %.15f \t %.15f \t %.15f \t %.15f\n", te, ye[0], ye[1], ye[2], ye[3]);

//             event_states.push_back(vector<double>());
//             event_states[l].push_back(te);
//             for ( int dimcnt = 0; dimcnt < dim; dimcnt++ ){
//                 event_states[l].push_back(ye[dimcnt]);
//             }
//             l += 1;
//         } 


//         trajectory[k][0] = t;
//         for (int dim_cnt = 1; dim_cnt < dim + 1; dim_cnt++ ){
//             trajectory[k][dim_cnt] = y[dim_cnt-1];
//         }

//         // printf("%.5e : \t", t);
//         // for (int dim_cnt = 0; dim_cnt < dim; dim_cnt++ ){
//         //     printf("%.5e \t", y[dim_cnt]);     
//         // }
//         // printf("\n");

//         // if ( k % 10000 == 0 )
//         //     printf(" Number of time steps: %5d\n", k);

            
//         fflush(stdout);

//     }
    
//     timer = clock() - timer;
//     printf ("It took me %lf secs. \n", ((double)timer)/CLOCKS_PER_SEC);

//     // Writing trajectory data to file
//     ofstream trajDataOut;
//     trajDataOut.open("Barbanis_traj_finalT"+to_string(int(tEnd))+"_E"+to_string(float(E))+".txt", ofstream::out | ofstream::trunc);
//     // trajDataOut.open("dm_traj_finalT"+to_string(int(tmax))+"_b"+to_string(int(mu))+".txt", ofstream::out );
//     for (int i=0; i< time_steps + 1; i++){
//         trajDataOut << std::fixed << std::setprecision(10) << setw(10) << trajectory[i][0] << "\t";
//         for (int k=1; k < dim + 1; k++)
//             trajDataOut << std::fixed << std::setprecision(15) << setw(15) << trajectory[i][k] << "\t";
//         trajDataOut << endl;
//     }
//     trajDataOut.close();


//     // Writing events to file
//     // printf("Size %d x %d\n ", l, event_states[0].size());
//     printf("Number of events : %d\n", l);

//     ofstream eventDataOut;
//     eventDataOut.open("event_finalT"+to_string(int(tEnd))+"_E"+to_string(float(E))+".txt", ofstream::out | ofstream::trunc);
//     for (int i=0; i< l; i++){
//         for (int k=0; k < dim + 1; k++)
//             eventDataOut << std::fixed << std::setprecision(15) << setw(15) << event_states[i][k] << "\t";
//         eventDataOut << endl;
//     }
//     eventDataOut.close();

//     return GSL_SUCCESS;
// }



// /*
// INTEGRATING USING VELOCITY VERLET ALGORITHM 
// */ 
// void vel_verlet_evolve(const double h, const vector<double> yi, vector<double>& yf){

//     double dV_dx[2]; // gradient of the potential function
//     double g1, g2, dg1_dx, dg2_dx, dg1_dy, dg2_dy;

//     // Calculate the force components at the current time 
//     dV_dx[0] = ( pow(OMEGA_X,2.0)*yi[0] + DELTA*pow(yi[1],2.0)  );
//     dV_dx[1] = ( pow(OMEGA_Y,2.0)*yi[1] + 2.0*DELTA*yi[0]*yi[1]  );

//     // Calculate position at new time
//     for ( int dimcnt = 0; dimcnt < dim/2; dimcnt++ ){
//         yf[dimcnt] = yi[dimcnt] + yi[dimcnt + 2]*h + (1/(2*M))*(-dV_dx[dimcnt])*pow(h,2);
//     }

//     // Calculate velocities at half-time step
//     for (int dimcnt = dim/2; dimcnt < dim; dimcnt++ ){
//         yf[dimcnt] = yi[dimcnt] + (1/(2*M))*(-dV_dx[dimcnt-2])*h;
//     }

//     // Calculate the force components at new time 
//     dV_dx[0] = ( pow(OMEGA_X,2.0)*yi[0] + DELTA*pow(yi[1],2.0)  );
//     dV_dx[1] = ( pow(OMEGA_Y,2.0)*yi[1] + 2.0*DELTA*yi[0]*yi[1]  ); 

//     // Calculate velocity at new time
//     for (int dimcnt = dim/2; dimcnt < dim; dimcnt++ ){
//         yf[dimcnt] = yf[dimcnt] + (1/(2*M))*(-dV_dx[dimcnt-2])*h;
//     }

// }

// /****************** EVENT LOCATION USING VELOCITY VERLET INTEGRATOR ****************************/ 
// void locate_event_vv(const double ti, const double tf, vector<double> yi, vector<double> yf, double& te, vector<double>& ye){

//     // int dim = 4;
//     double eventTime;
//     double h = 1e-6;
//     // double h;

//     // Bisection method parameters
//     int iter = 0;
//     int iterMax = 100;

//     // printf("Event is between %.5f \t %.5f \n", ti, tf);
//     // for ( int dim_cnt = 0; dim_cnt < dim; dim_cnt++ ){
//     //     printf("%.5f \t %.5f \n", yi[dim_cnt], yf[dim_cnt]);    
//     // }
    
//     double prev_val = eventfun(yi);
//     double curr_val = eventfun(yf);

//     // cout << prev_val << " " << curr_val << endl;

//     double tLeft = ti;
//     double tRight = tf;
//     double ev_t_guess, tInitial;
//     vector<double> y(dim);

//     double val_mid;
//     // double params = B;

//     vector<double> y_start = yi; 
//     while ( ( iter < iterMax )  ||  ( fabs(prev_val - curr_val) > 1e-10 ) ){

//         ev_t_guess = 0.5*( tLeft + tRight );
//         // printf("New guess for event time : %.10f\n", ev_t_guess);
//         tInitial = tLeft;
//         h = ev_t_guess - tLeft;

//         // double *y_start = &yMid[0]; // because yMid is STL vector
        

//         // printf("%f \t %f \n", tInitial, ev_t_guess);
//         // printf("States at %f ... \n ", tInitial);
//         // for ( int dimcnt=0; dimcnt<dim;dimcnt++){
//         //     printf("%f \t", y_start[dimcnt]);
//         // }
//         // printf("\n");


//         vel_verlet_evolve(h, y_start, y);

//         // printf("States at %f ... \n ", tInitial + h);
//         // for ( int dimcnt=0; dimcnt<dim;dimcnt++){
//         //     printf("%f \t", y_start[dimcnt]);
//         // }
//         // printf("\n");

//         // event value at the mid point
//         val_mid = eventfun(y); 
        
//         // printf("Event value is %.15f at guess %.15f \n", val_mid, ev_t_guess);

//         if ( val_mid/fabs(val_mid) == curr_val/fabs(curr_val) ){ // move the curr point to mid point

//             tRight = ev_t_guess;
//             curr_val = val_mid;
//             yf = y;
//             // y_start = yi;
//             // printf("Event is between %f and %f\n", tLeft, tRight);

//         } 
//         else if ( val_mid/fabs(val_mid) == prev_val/fabs(prev_val) ){ // move the prev point to mid point
            
//             tLeft = ev_t_guess;
//             prev_val = val_mid;
//             y_start = y;
//             // printf("Event is between %f and %f .\n", tLeft, tRight);
            
//         }

//        eventTime = ev_t_guess;

//         iter += 1;
//     }

//     te = eventTime;
//     ye = y;
//     // printf("%f \t %f \t %f \t %f\n", ye[0], ye[1], ye[2], ye[3]);
    
// }

