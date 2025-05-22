/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   March 5, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rdt dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
  v[RHO] = 1.0;
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  #if HAVE_ENERGY
  v[PRS] = 1.0;
  #endif
  v[TRC] = 0.0;

  #if PHYSICS == MHD || PHYSICS == RMHD
  v[BX1] = 0.0;
  v[BX2] = 0.0;
  v[BX3] = 0.0;

  v[AX1] = 0.0;
  v[AX2] = 0.0;
  v[AX3] = 0.0;
  #endif
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
  int i, j, k, nv;
  g_gamma = 5/3.;
  g_minCoolingTemp = 1.0e+04;
  
  TOT_LOOP(k, j, i) {
    d->Vc[RHO][k][j][i] = 1.0;
    d->Vc[PRS][k][j][i] = ((d->Vc[RHO][k][j][i]*UNIT_DENSITY)/(0.609*CONST_mp) * CONST_kB * g_inputParam[TINI])/(UNIT_DENSITY*pow(UNIT_VELOCITY,2));
    d->Vc[VX1][k][j][i] = 1.;
    d->Vc[VX2][k][j][i] = 0.;
    d->Vc[VX3][k][j][i] = 0.;
    #if COOLING==GRACKLE
    double tiny_number = 1.e-20;
    double XH_mass = 0.715768377353088514;
    double ZmetSol_mass = 0.01295; 
    d->Vc[TEMP][k][j][i] = g_inputParam[TINI];
    d->Vc[MU][k][j][i] = 0.609;
    NIONS_LOOP(nv) d->Vc[nv][k][j][i] = 0;
    d->Vc[X_HI][k][j][i]     = tiny_number;
    d->Vc[X_HII][k][j][i]    = XH_mass;
    d->Vc[Y_HeI][k][j][i]    = tiny_number;
    d->Vc[Y_HeII][k][j][i]   = tiny_number;
    d->Vc[Y_HeIII][k][j][i]  = (1-XH_mass);
    d->Vc[X_HM][k][j][i]     = tiny_number;
    d->Vc[X_H2I][k][j][i]    = tiny_number;
    d->Vc[X_H2II][k][j][i]   = tiny_number;
    d->Vc[X_DI][k][j][i]     = tiny_number;
    d->Vc[X_DII][k][j][i]    = 2.0 * 3.4e-05;
    d->Vc[X_HDI][k][j][i]    = tiny_number;
    d->Vc[elec][k][j][i]     = (d->Vc[X_HII][k][j][i] + d->Vc[X_DII][k][j][i] + 
	                           (d->Vc[Y_HeII][k][j][i]+2*d->Vc[Y_HeIII][k][j][i])/4.)*d->Vc[RHO][k][j][i];
    d->Vc[Z_MET][k][j][i]    = 1.0;
    #endif
  }
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{
  int k, j, i;
  
  double temp = 0., mass = 0., energy = 0., vol = 0.;
  static double energy_old = 0.;
  static double temperature_old = 0.;
  static double time_old = 0.;
  if (g_stepNumber==0) {
    temperature_old = g_inputParam[TINI];
    energy_old = (d->Vc[PRS][0][0][0]/d->Vc[RHO][0][0][0])/(g_gamma-1);
  }
  DOM_LOOP(k, j, i) {
    temp   += (d->Vc[RHO][k][j][i]*d->Vc[TEMP][k][j][i]*grid->dV[k][j][i]);
    mass   += (d->Vc[RHO][k][j][i]*grid->dV[k][j][i]);
    energy += (d->Vc[PRS][k][j][i]*grid->dV[k][j][i]/(g_gamma-1));
    vol    += grid->dV[k][j][i];
  }
  #ifdef PARALLEL
  int transfer_size = 4;
  int transfer = 0;
  double sendArray[transfer_size], recvArray[transfer_size];
  sendArray[transfer++] = temp; sendArray[transfer++] = mass; sendArray[transfer++] = energy; sendArray[transfer++] = vol;
  MPI_Allreduce (sendArray, recvArray, transfer_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  transfer = 0;
  temp = recvArray[transfer++]; mass = recvArray[transfer++]; energy = recvArray[transfer++]; vol = recvArray[transfer++];
  #endif
  double T_avg = temp/mass;
  double energy_avg = energy/vol;
  double dE = 0., dT = 0.;
  double lambda = 0.;
  if (g_stepNumber>0) {
    dE = energy_avg - energy_old;
    dT = T_avg - temperature_old;
    double derv = dE/(g_time-time_old);
    derv = derv * UNIT_DENSITY * pow(UNIT_VELOCITY, 3) * pow(UNIT_LENGTH, -1);
    double ndens = (mass/vol)*UNIT_DENSITY/(0.609*CONST_mp);
    double nH = ndens*0.609*0.71;
    lambda = -derv/pow(nH, 2);
    time_old = g_time;
    energy_old = energy_avg;
    temperature_old = T_avg;
  }
  /* ---- Write ascii file "analysis.dat" to disk ---- */
  if (prank == 0){
    char fname[512];
    FILE *fp;
    sprintf (fname, "%s/analysis.dat",RuntimeGet()->output_dir);
    if (g_stepNumber == 0) { /* Open for writing only when we are starting */
      fp = fopen(fname,"w"); /* from beginning */
      fprintf (fp,"# %s\t=\t%.5e\n", "UNIT_LENGTH", UNIT_LENGTH);
      fprintf (fp,"# %s\t=\t%.5e\n", "UNIT_DENSITY", UNIT_DENSITY);
      fprintf (fp,"# %s\t=\t%.5e\n", "UNIT_VELOCITY", UNIT_VELOCITY);
      // Header
      fprintf (fp,"# (1)%s\t(2)%s\t(3)%s\t\(4)%s\t\n",
               "time (code)", "temperature (K)", "energy (erg/cm^3)", "lambda (ergcm^3/s)");
      fclose(fp);
    }
    /* Append numeric data */
    fp = fopen(fname,"a");
    fprintf (fp, "%12.6e\t\t%12.6e\t\t%12.6e\t\t%12.6e\t\t\n",
             g_time, T_avg, (energy_avg*UNIT_DENSITY*pow(UNIT_VELOCITY,2)*pow(UNIT_LENGTH,-3)), lambda);
    fclose(fp);
  }
  /* --- end of writing "analyis.dat" file --- */
}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){}
  }

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X1_END){  /* -- X1_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_END){  /* -- X2_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_END){  /* -- X3_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif
