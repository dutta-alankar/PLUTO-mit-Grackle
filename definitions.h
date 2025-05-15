#define  PHYSICS                        HD
#define  DIMENSIONS                     1
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        GRACKLE
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            1

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  TINI                           0

/* [Beg] user-defined constants (do not change this line) */
#define  UNIT_DENSITY                   (1.0e-03*0.609*CONST_mp)
#define  UNIT_LENGTH                    CONST_pc
#define  UNIT_VELOCITY                  1.0e+05

/* [End] user-defined constants (do not change this line) */

// #define SHOW_TIMING                     YES
#define SHOW_TIME_STEPS                 YES

