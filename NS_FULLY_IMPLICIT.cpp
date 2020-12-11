/* Notation for this code

u           =   velocity component in x direction for current Time Step
u_star      =   Predicted velocity component in x direction
u_Next      =   velocity component in x direction for Next Time Step
v           =   velocity component in y direction for current Time Step
v_star      =   Predicted velocity component in y direction
v_Next      =   velocity component in y direction for Next Time Step
p_Next      =   vorticity for next time step                                */

#include <iostream>
#include<math.h>
#include <fstream>

using namespace std;

int i,j;
const int Nx = 50, Ny = 50;
double Re = 100, L = 1,dt = 0.001;
double dx,dy,Vp;

void Gauss_Seidal_Predictor(double u[Nx+2][Ny+2],double u_star[Nx+2][Ny+2],double v[Nx+2][Ny+2],double v_star[Nx+2][Ny+2],double p_Next[Nx+2][Ny+2]);   // Predictor function
void Gauss_Seidal_uStar(double u[Nx+2][Ny+2],double u_star[Nx+2][Ny+2],double v[Nx+2][Ny+2],double v_star[Nx+2][Ny+2],double p_star[Nx+2][Ny+2]);
void Gauss_Seidal_Pressure(double u_star[Nx+2][Ny+2],double v_star[Nx+2][Ny+2],double p_Next[Nx+2][Ny+2]);                                              // Gauss Seidal function for Pressure
void Gauss_Seidal_pStar(double u_star[Nx+2][Ny+2],double v_star[Nx+2][Ny+2],double p_star[Nx+2][Ny+2]);
void Gauss_Seidal_Corrector( double u[Nx+2][Ny+2],double u_Next[Nx+2][Ny+2],double u_star[Nx+2][Ny+2],double v[Nx+2][Ny+2],double v_Next[Nx+2][Ny+2],double v_star[Nx+2][Ny+2],double p_Next[Nx+2][Ny+2]);   // Corrector function
void exporting_files(double u[Nx+2][Ny+2],double v[Nx+2][Ny+2],double p_Next[Nx+2][Ny+2]);                                                              // function for exporting data

int main()
{
    int z = 0;

    double u[Nx+2][Ny+2] = {0}, u_star[Nx+2][Ny+2] = {0} , u_Next[Nx+2][Ny+2] = {0};
    double v[Nx+2][Ny+2] = {0}, v_star[Nx+2][Ny+2] = {0} , v_Next[Nx+2][Ny+2] = {0};
    double p_star[Nx+2][Ny+2] = {0} ,p_Next[Nx+2][Ny+2] = {0} , p_star_old[Nx+2][Ny+2] = {0} ;

    double res,err;

    double error = 1;

  while ( error > 1e-03 )
  {
        z = z+1;
        cout<<z<<endl;

        err = 1;

        while ( err > 1e-06 )
        {


            Gauss_Seidal_uStar( u,u_star,v,v_star,p_star );                         // Function for Predictor Step

            Gauss_Seidal_pStar( u_star,v_star,p_star );                             // Function for Gauss Seidal Pressure

            res = 0;

            for ( j = 0 ; j <= Ny+1 ; j++ )
            {
                for ( i = 0 ; i <= Nx+1 ; i++ )
                {
                    res =  res + pow((p_star[i][j] - p_star_old[i][j]),2);
                    p_star_old[i][j] = p_star[i][j];
                }
            }
                    err = sqrt(res/(Nx*Ny));

          }

        Gauss_Seidal_Corrector( u,u_Next,u_star,v,v_Next,v_star,p_Next );            // Function for Corerctor Step


    // Checking for Steady State

    error = 0;

    for ( j = 0 ; j <= Ny+1 ; j++ )
    {
        for ( i = 0 ; i <= Nx+1 ; i++ )
        {
            error = error + pow((u_Next[i][j] - u[i][j]),2);            // Residue for difference in current and next time step values
            u[i][j]  = u_Next[i][j];                                    // updating values for next time step
            v[i][j]  = v_Next[i][j];
        }
    }

    error = sqrt(error/(Nx*Ny));
    error = error/dt;
  }

    exporting_files(u,v,p_Next);

    return 0;
}

// Function for ustar_star

void Gauss_Seidal_uStar(double u[Nx+2][Ny+2],double u_star[Nx+2][Ny+2],double v[Nx+2][Ny+2],double v_star[Nx+2][Ny+2],double p_star[Nx+2][Ny+2])
{
    dx = L/Nx;
    dy = L/Ny;
    Vp = dx*dy;

    double fe[Nx+2][Ny+2], fw[Nx+2][Ny+2], fn[Nx+2][Ny+2], fs[Nx+2][Ny+2], uff, vff, fduf, fdvf;
    double ue, uw, un, us;
    double ve, vw, vn, vs;

    for ( j = 1 ; j < Ny+1 ; j++ )
    {
        for ( i = 1 ; i < Nx+1 ; i++ )
        {
            fe[i][j] =   0.5 * ( u_star[i][j] + u_star[i+1][j] ) * dy  -  dt * ( p_star[i+1][j] - p_star[i][j]   ) * dy/dx;              // Calculating volume flux on all faces
            fw[i][j] = - 0.5 * ( u_star[i][j] + u_star[i-1][j] ) * dy  +  dt * ( p_star[i][j]   - p_star[i-1][j] ) * dy/dx;
            fn[i][j] =   0.5 * ( v_star[i][j] + v_star[i][j+1] ) * dx  -  dt * ( p_star[i][j+1] - p_star[i][j]   ) * dx/dy;
            fs[i][j] = - 0.5 * ( v_star[i][j] + v_star[i][j-1] ) * dx  +  dt * ( p_star[i][j]   - p_star[i][j-1] ) * dx/dy;
        }
    }

    double err,x;
    double res1,res2,r1,r2;

            err = 1;

    while ( err > 1e-06 )
    {
            res1 = 0;
            res2 = 0;

        for ( j = 1 ; j < Ny+1 ; j++ )
        {
            for ( i = 1 ; i < Nx+1 ; i++ )
            {
                        x = 0;
                if ( fe[i][j] >= 0) {         // Upwind scheme for Convected term for boundary cells since we cant use QUICK at boundary
                        ue = u_star[i][j];
                        ve = v_star[i][j];
                        x = x + fe[i][j];}
                else{
                        ue = u_star[i+1][j];
                        ve = v_star[i+1][j];}

                if ( fw[i][j] >= 0){
                        uw = u_star[i][j];
                        vw = v_star[i][j];
                        x = x + fw[i][j];}
                else{
                        uw = u_star[i-1][j];
                        vw = v_star[i-1][j];}

                if ( fn[i][j] >= 0){
                        un = u_star[i][j];
                        vn = v_star[i][j];
                        x = x + fn[i][j];}
                else{
                        un = u_star[i][j+1];
                        vn = v_star[i][j+1];}

                if ( fs[i][j] >= 0){
                        us = u_star[i][j];
                        vs = v_star[i][j];
                        x = x + fs[i][j];}
                else{
                        us = u_star[i][j-1];
                        vs = v_star[i][j-1];}

                        if ( i > 1 && j > 1 )                                            // Applying QUICK scheme for conveted terms for interior cells
                        {
                            if ( i < Nx && j < Ny )
                            {
                                if ( fe[i][j] >= 0){
                                    ue = 6.0/8*u_star[i][j] + 3.0/8*u_star[i+1][j] - 1.0/8*u_star[i-1][j];
                                    ve = 6.0/8*v_star[i][j] + 3.0/8*v_star[i+1][j] - 1.0/8*v_star[i-1][j];}
                                else{
                                    ue = 6.0/8*u_star[i+1][j] + 3.0/8*u_star[i][j] - 1.0/8*u_star[i+2][j];
                                    ve = 6.0/8*v_star[i+1][j] + 3.0/8*v_star[i][j] - 1.0/8*v_star[i+2][j];}

                                if ( fw[i][j] >= 0){
                                    uw = 6.0/8*u_star[i-1][j] + 3.0/8*u_star[i][j] - 1.0/8*u_star[i-2][j];
                                    vw = 6.0/8*v_star[i-1][j] + 3.0/8*v_star[i][j] - 1.0/8*v_star[i-2][j];}
                                else{
                                    uw = 6.0/8*u_star[i][j] + 3.0/8*u_star[i-1][j] - 1.0/8*u_star[i+1][j];
                                    vw = 6.0/8*v_star[i][j] + 3.0/8*v_star[i-1][j] - 1.0/8*v_star[i+1][j];}

                                if ( fn[i][j] >= 0){
                                    un = 6.0/8*u_star[i][j] + 3.0/8*u_star[i][j+1] - 1.0/8*u_star[i][j-1];
                                    vn = 6.0/8*v_star[i][j] + 3.0/8*v_star[i][j+1] - 1.0/8*v_star[i][j-1];}
                                else{
                                    un = 6.0/8*u_star[i][j+1] + 3.0/8*u_star[i][j] - 1.0/8*u_star[i][j+2];
                                    vn = 6.0/8*v_star[i][j+1] + 3.0/8*v_star[i][j] - 1.0/8*v_star[i][j+2];}

                                if ( fs[i][j] >= 0){
                                    us = 6.0/8*u_star[i][j-1] + 3.0/8*u_star[i][j] - 1.0/8*u_star[i][j-2];
                                    vs = 6.0/8*v_star[i][j-1] + 3.0/8*v_star[i][j] - 1.0/8*v_star[i][j-2];}
                                else{
                                    us = 6.0/8*u_star[i][j] + 3.0/8*u_star[i][j-1] - 1.0/8*u_star[i][j+1];
                                    vs = 6.0/8*v_star[i][j] + 3.0/8*v_star[i][j-1] - 1.0/8*v_star[i][j+1];}
                            }
                        }

            uff = fe[i][j] * ue + fw[i][j] * uw + fn[i][j] * un + fs[i][j] * us ;                        // Calculating total volume Flux
            vff = fe[i][j] * ve + fw[i][j] * vw + fn[i][j] * vn + fs[i][j] * vs ;

            fduf  = ( u_star[i+1][j] -2*u_star[i][j] + u_star[i-1][j] )/( dx * dx) + ( u_star[i][j+1] -2*u_star[i][j] + u_star[i][j-1] )/( dy * dy) ;
            fdvf  = ( v_star[i+1][j] -2*v_star[i][j] + v_star[i-1][j] )/( dx * dx) + ( v_star[i][j+1] -2*v_star[i][j] + v_star[i][j-1] )/( dy * dy) ;

            r1 = Vp/dt * ( u[i][j] - u_star[i][j] ) - uff + Vp/Re * fduf ;
            res1 = res1 + r1 * r1;
            u_star[i][j] = r1/( Vp/dt + x + Vp/Re*( 2/(dx*dx) + 2/(dy*dy) ) ) + u_star[i][j] ;

            r2 = Vp/dt * ( v[i][j] - v_star[i][j] ) - vff + Vp/Re * fdvf ;
            res2 = res2 + r2 *r2;
            v_star[i][j] = r2/( Vp/dt + x + Vp/Re*( 2/(dx*dx) + 2/(dy*dy) ) ) + v_star[i][j] ;

            }
        }
            err = sqrt( (res1+res2)/(Nx*Ny));
    }



    // BC for u_star velocities

    for ( j = 1 ; j <= Ny ; j++ )               // Left
    {
            u_star[0][j] = - u_star[1][j];
            v_star[0][j] = - v_star[1][j];

    }
    for ( j = 1 ; j <= Ny ; j++ )               // Right
    {
            u_star[Nx+1][j] = - u_star[Nx][j];
            v_star[Nx+1][j] = - v_star[Nx][j];

    }
    for ( i = 1 ; i <= Nx ; i++ )               // Bottom
    {
            u_star[i][0] = - u_star[i][1];
            v_star[i][0] = - v_star[i][1];

    }
    for ( i = 1 ; i <= Nx ; i++ )               // Top
    {
            u_star[i][Ny+1] = 2 - u_star[i][Ny];
            v_star[i][Ny+1] =   - v_star[i][Ny];
    }
}

// Gauss Seidal Function for pStar

void Gauss_Seidal_pStar(double u_star[Nx+2][Ny+2],double v_star[Nx+2][Ny+2],double p_star[Nx+2][Ny+2])
{
    dx = L/Nx;
    dy = L/Ny;
    Vp = dx*dy;

    double R,er,res;

    er  = 1;

    while ( er > 1e-06 )                                    // error loop for gauss seidal stream function
    {
        res = 0;

        for ( j = 1 ; j < Ny+1 ; j++ )
        {
            for ( i = 1 ; i < Nx+1 ; i++ )
            {
                R  = 0.5 * ( (u_star[i+1][j] - u_star[i-1][j])/dx + (v_star[i][j+1] - v_star[i][j-1])/dy ) - dt * ( p_star[i+1][j] - 2*p_star[i][j] + p_star[i-1][j])/(dx*dx) - dt * ( p_star[i][j+1] - 2*p_star[i][j] + p_star[i][j-1])/(dy*dy) ;
                res = res + R * R;                                                     // Adding square of residues for RMS
                p_star[i][j] = - 1.5*R/( 2*dt/(dx*dx) + 2*dt/(dy*dy) ) + p_star[i][j];           // correction for next iteration
            }
        }
        // BC for p_Next

        for ( j = 1 ; j <= Ny ; j++ )               // Left
        {
            p_star[0][j] = p_star[1][j];
        }
        for ( j = 1 ; j <= Ny ; j++ )               // Right
        {
            p_star[Nx+1][j] = p_star[Nx][j];
        }
        for ( i = 1 ; i <= Nx ; i++ )               // Bottom
        {
            p_star[i][0] = p_star[i][1];
        }
        for ( i = 1 ; i <= Nx ; i++ )               // Top
        {
            p_star[i][Ny+1] = p_star[i][Ny];
        }

            er = sqrt(res/(Nx*Ny));
    }

}


// Function for Corrector Step

void Gauss_Seidal_Corrector( double u[Nx+2][Ny+2],double u_Next[Nx+2][Ny+2],double u_star[Nx+2][Ny+2],double v[Nx+2][Ny+2],double v_Next[Nx+2][Ny+2],double v_star[Nx+2][Ny+2],double p_Next[Nx+2][Ny+2])
{
    dx = L/Nx;
    dy = L/Ny;
    Vp = dx*dy;

    double fe, fw, fn, fs, uff, vff, fduf, fdvf;
    double ue, uw, un, us;
    double ve, vw, vn, vs;

    double err,x;
    double res1,res2,r1,r2;

            err = 1;

    while ( err > 1e-06 )
    {
            res1 = 0;
            res2 = 0;

        for ( j = 1 ; j < Ny+1 ; j++ )
        {
            for ( i = 1 ; i < Nx+1 ; i++ )
            {
                x = 0;
                fe =   0.5 * ( u_star[i][j] + u_star[i+1][j] ) * dy  -  dt * ( p_Next[i+1][j] - p_Next[i][j]   ) * dy/dx;              // Calculating volume flux on all faces
                fw = - 0.5 * ( u_star[i][j] + u_star[i-1][j] ) * dy  +  dt * ( p_Next[i][j]   - p_Next[i-1][j] ) * dy/dx;
                fn =   0.5 * ( v_star[i][j] + v_star[i][j+1] ) * dx  -  dt * ( p_Next[i][j+1] - p_Next[i][j]   ) * dx/dy;
                fs = - 0.5 * ( v_star[i][j] + v_star[i][j-1] ) * dx  +  dt * ( p_Next[i][j]   - p_Next[i][j-1] ) * dx/dy;

                if ( fe >= 0) {         // Upwind scheme for Convected term for boundary cells since we cant use QUICK at boundary
                        ue = u_Next[i][j];
                        ve = v_Next[i][j];
                        x = x + fe;}
                else{
                        ue = u_Next[i+1][j];
                        ve = v_Next[i+1][j];}

                if ( fw >= 0){
                        uw = u_Next[i][j];
                        vw = v_Next[i][j];
                        x = x + fw;}
                else{
                        uw = u_Next[i-1][j];
                        vw = v_Next[i-1][j];}

                if ( fn >= 0){
                        un = u_Next[i][j];
                        vn = v_Next[i][j];
                        x = x + fn;}
                else{
                        un = u_Next[i][j+1];
                        vn = v_Next[i][j+1];}

                if ( fs >= 0){
                        us = u_Next[i][j];
                        vs = v_Next[i][j];
                        x = x + fs;}
                else{
                        us = u_Next[i][j-1];
                        vs = v_Next[i][j-1];}

                        if ( i > 1 && j > 1 )                                            // Applying QUICK scheme for conveted terms for interior cells
                        {
                            if ( i < Nx && j < Ny )
                            {
                                if ( fe >= 0){
                                    ue = 6.0/8*u_Next[i][j] + 3.0/8*u_Next[i+1][j] - 1.0/8*u_Next[i-1][j];
                                    ve = 6.0/8*v_Next[i][j] + 3.0/8*v_Next[i+1][j] - 1.0/8*v_Next[i-1][j];}
                                else{
                                    ue = 6.0/8*u_Next[i+1][j] + 3.0/8*u_Next[i][j] - 1.0/8*u_Next[i+2][j];
                                    ve = 6.0/8*v_Next[i+1][j] + 3.0/8*v_Next[i][j] - 1.0/8*v_Next[i+2][j];}

                                if ( fw >= 0){
                                    uw = 6.0/8*u_Next[i-1][j] + 3.0/8*u_Next[i][j] - 1.0/8*u_Next[i-2][j];
                                    vw = 6.0/8*v_Next[i-1][j] + 3.0/8*v_Next[i][j] - 1.0/8*v_Next[i-2][j];}
                                else{
                                    uw = 6.0/8*u_Next[i][j] + 3.0/8*u_Next[i-1][j] - 1.0/8*u_Next[i+1][j];
                                    vw = 6.0/8*v_Next[i][j] + 3.0/8*v_Next[i-1][j] - 1.0/8*v_Next[i+1][j];}

                                if ( fn >= 0){
                                    un = 6.0/8*u_Next[i][j] + 3.0/8*u_Next[i][j+1] - 1.0/8*u_Next[i][j-1];
                                    vn = 6.0/8*v_Next[i][j] + 3.0/8*v_Next[i][j+1] - 1.0/8*v_Next[i][j-1];}
                                else{
                                    un = 6.0/8*u_Next[i][j+1] + 3.0/8*u_Next[i][j] - 1.0/8*u_Next[i][j+2];
                                    vn = 6.0/8*v_Next[i][j+1] + 3.0/8*v_Next[i][j] - 1.0/8*v_Next[i][j+2];}

                                if ( fs >= 0){
                                    us = 6.0/8*u_Next[i][j-1] + 3.0/8*u_Next[i][j] - 1.0/8*u_Next[i][j-2];
                                    vs = 6.0/8*v_Next[i][j-1] + 3.0/8*v_Next[i][j] - 1.0/8*v_Next[i][j-2];}
                                else{
                                    us = 6.0/8*u_Next[i][j] + 3.0/8*u_Next[i][j-1] - 1.0/8*u_Next[i][j+1];
                                    vs = 6.0/8*v_Next[i][j] + 3.0/8*v_Next[i][j-1] - 1.0/8*v_Next[i][j+1];}
                            }
                        }

            uff = fe * ue + fw * uw + fn * un + fs * us ;                        // Calculating total volume Flux
            vff = fe * ve + fw * vw + fn * vn + fs * vs ;

            fduf  = ( u_Next[i+1][j] -2*u_Next[i][j] + u_Next[i-1][j] )/( dx * dx) + ( u_Next[i][j+1] -2*u_Next[i][j] + u_Next[i][j-1] )/( dy * dy) ;
            fdvf  = ( v_Next[i+1][j] -2*v_Next[i][j] + v_Next[i-1][j] )/( dx * dx) + ( v_Next[i][j+1] -2*v_Next[i][j] + v_Next[i][j-1] )/( dy * dy) ;

            r1 = u[i][j] - u_Next[i][j] - dt/Vp * uff + dt/Re * fduf - ( dt * 0.5/dx * ( p_Next[i+1][j] - p_Next[i-1][j] ) );
            res1 = res1 + r1 * r1;
            u_Next[i][j] = r1/( 1 + dt/Vp * x + dt/Re*( 2/(dx*dx) + 2/(dy*dy) ) ) + u_Next[i][j] ;

            r2 = v[i][j] - v_Next[i][j] - dt/Vp * vff + dt/Re * fdvf - ( dt * 0.5/dy * ( p_Next[i][j+1] - p_Next[i][j-1] ) );
            res2 = res2 + r2 *r2;
            v_Next[i][j] = r2/( 1 + dt/Vp * x + dt/Re*( 2/(dx*dx) + 2/(dy*dy) ) ) + v_Next[i][j] ;

            }
        }
            err = sqrt( (res1+res2)/(Nx*Ny));
    }


    // BC for u_Next velocities
    for ( j = 1 ; j <= Ny ; j++ )               // Left
    {
            u_Next[0][j] = - u_Next[1][j];
            v_Next[0][j] = - v_Next[1][j];

    }
    for ( j = 1 ; j <= Ny ; j++ )               // Right
    {
            u_Next[Nx+1][j] = - u_Next[Nx][j];
            v_Next[Nx+1][j] = - v_Next[Nx][j];

    }
    for ( i = 1 ; i <= Nx ; i++ )               // Bottom
    {
            u_Next[i][0] = - u_Next[i][1];
            v_Next[i][0] = - v_Next[i][1];

    }
    for ( i = 1 ; i <= Nx ; i++ )               // Top
    {
            u_Next[i][Ny+1] = 2 - u_Next[i][Ny];
            v_Next[i][Ny+1] =   - v_Next[i][Ny];
    }

}

// Exporting Data in .dat Files

void exporting_files(double u[Nx+2][Ny+2],double v[Nx+2][Ny+2],double p_Next[Nx+2][Ny+2])
{
    // Exporting the Velocities

    double U[52],V[52];

    ofstream myfile;
    myfile.open ("U_NS_QUICK_FULL_IMPLICIT.dat");           // Center Line U velocity_calculation

        U[0] = 0.5 * ( ( u[26][1] + u[25][1] )/2  + ( u[26][0] + u[25][0] )/2 );
        myfile <<U[0]<<endl;
    for ( i = 1 ; i <= Ny ; i++ )
    {
        U[i] = ( u[26][i] + u[25][i] )/2  ;
        myfile <<U[i]<<endl;
    }
        U[Ny+1] = 0.5 * ( ( u[26][Ny+1] + u[25][Ny+1] )/2  + ( u[26][Ny] + u[25][Ny] )/2 );
        myfile <<U[Ny+1]<<endl;
    myfile.close();

    myfile.open ("V_NS_QUICK_FULL_IMPLICIT.dat");           // Center line V velocity_calculation

        V[0] = 0.5 * ( ( v[1][26] + v[1][25] )/2 + ( v[0][26] + v[0][25] )/2 );
        myfile << V[0]<<endl;
    for ( i = 1 ; i <= Nx ; i++ )
    {
        V[i] = ( v[i][26] + v[i][25] )/2 ;
        myfile << V[i]<<endl;
    }
        V[Nx+1] = 0.5 * ( ( v[Nx+1][26] + v[Nx+1][25] )/2 + ( v[Nx][26] + v[Nx][25] )/2 );
        myfile << V[Nx+1]<<endl;
    myfile.close();

    for ( i = 1 ; i <= Nx ; i++ )
    {
        u[i][0] = ( u[i][1] + u[i][0] )/2 ;                         //Bottom
        v[i][0] = ( v[i][1] + v[i][0] )/2 ;
        p_Next[i][0] = ( p_Next[i][1] + p_Next[i][0] )/2 ;

        u[i][Ny+1] = ( u[i][Ny+1] + u[i][Ny] )/2 ;                  // Top
        v[i][Ny+1] = ( v[i][Ny+1] + v[i][Ny] )/2 ;
        p_Next[i][Ny+1] = ( p_Next[i][Ny+1] + p_Next[i][Ny] )/2 ;
    }

    for ( j = 1 ; i <= Ny ; j++ )
    {
        u[0][j] = ( u[1][j] + u[0][j] )/2 ;                         // Left
        v[0][j] = ( v[1][j] + v[0][j] )/2 ;
        p_Next[0][j] = ( p_Next[1][j] + p_Next[0][j] )/2 ;

        u[Nx+1][j] = ( u[Nx+1][j] + u[Nx][j] )/2 ;                 // Right
        v[Nx+1][j] = ( v[Nx+1][j] + v[Nx][j] )/2 ;
        p_Next[Nx+1][j] = ( p_Next[Nx+1][j] + p_Next[Nx][j] )/2 ;
    }

        u[0][0] = 0.5 * ( ( u[1][1] + u[1][0] )/2 + ( u[0][0] + u[0][1] )/2 );                                  // BottomLeft corner
        v[0][0] = 0.5 * ( ( v[1][1] + v[1][0] )/2 + ( v[0][0] + v[0][1] )/2 );
        p_Next[0][0] = 0.5 * ( ( p_Next[1][1] + p_Next[1][0] )/2 + ( p_Next[0][0] + p_Next[0][1] )/2 );

        u[Nx+1][0] = 0.5 * ( ( u[Nx+1][1] + u[Nx+1][0] )/2 + ( u[Nx][0] + u[Nx][1] )/2 );                       // Bottom Right corner
        v[Nx+1][0] = 0.5 * ( ( v[Nx+1][1] + v[Nx+1][0] )/2 + ( v[Nx][0] + v[Nx][1] )/2 );
        p_Next[Nx+1][0] = 0.5 * ( ( p_Next[Nx+1][1] + p_Next[Nx+1][0] )/2 + ( p_Next[Nx][0] + p_Next[Nx][1] )/2 );

        u[0][Ny+1] = 0.5 * ( ( u[1][Ny+1] + u[1][Ny] )/2 + ( u[0][Ny+1] + u[0][Ny] )/2 );                       // Top Left corner
        v[0][Ny+1] = 0.5 * ( ( v[1][Ny+1] + v[1][Ny] )/2 + ( v[0][Ny+1] + v[0][Ny] )/2 );
        p_Next[0][Ny+1] = 0.5 * ( ( p_Next[1][Ny+1] + p_Next[1][Ny] )/2 + ( p_Next[0][Ny+1] + p_Next[0][Ny] )/2 );

        u[Nx+1][Ny+1] = 0.5 * ( ( u[Nx+1][Ny+1] + u[Nx+1][Ny] )/2 + ( u[Nx][Ny+1] + u[Nx][Ny] )/2 );            // Top Right corner
        v[Nx+1][Ny+1] = 0.5 * ( ( v[Nx+1][Ny+1] + v[Nx+1][Ny] )/2 + ( v[Nx][Ny+1] + v[Nx][Ny] )/2 );
        p_Next[Nx+1][Ny+1] = 0.5 * ( ( p_Next[Nx+1][Ny+1] + p_Next[Nx+1][Ny] )/2 + ( p_Next[Nx][Ny+1] + p_Next[Nx][Ny] )/2 );

        // Saving contour data

    myfile.open ("DATA_NS_QUICK_FULL_IMPLICIT.dat");

    for ( j = 0 ; j <= Nx+1 ; j++ ){
    for ( i = 0 ; i <= Nx+1 ; i++ )
    {
        myfile <<u[i][j]<<","<< v[i][j]<<","<< p_Next[i][j]<<endl;                                    // Saving Data
    }}
    myfile.close();

}
