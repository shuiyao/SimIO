#define MAXDIM 3
#define NMETALS 4
#define forever for(;;)

typedef float Real;

struct gas_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real rho;
    Real temp;
    Real hsmooth;
    Real metals ;
    Real phi ;
} ;

struct dark_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real eps;
    Real phi ;
} ;

struct star_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real metals ;
    Real tform ;
    Real eps;
    Real phi ;
} ;

struct dump {
    double time ;
    int nbodies ;
    int ndim ;
    int nsph ;
    int ndark ;
    int nstar ;
    int pad;
} ;

struct dump tipsy_header ;

struct aux_gas_data
{
  float metal[NMETALS];
  float sfr;
  float tmax;
  float delaytime;
  float ne;
  float nh;
  int nspawn;

#if !defined(JMG_DEBUG)
/* SHUIYAO: 13-04-09 */
  short int nrec;
#endif // JMG_DEBUG

#ifdef OUTPUT_ALPHA
  float alpha;
#endif
};

struct aux_star_data
{
  float metal[NMETALS];
  float age;
  float tmax;
  int nspawn;

#if !defined(JMG_DEBUG)
/* SHUIYAO: 13-04-09 */
  short int nrec;
#endif // JMG_DEBUG

};

#ifdef OUTPUT_TIPSY_AW
struct aw_gas_data
{
  short int wind_flag;
  float mass_cloud;
  float rho;
  float temp;
  float dtcool;
  float dudt;
  float rcloud;
#ifdef WIND_NGB_STAT
  int numngb_as_wind;
  int numngb_nonwind;
  float r_nearest_ngb;
#endif  
};
#endif
