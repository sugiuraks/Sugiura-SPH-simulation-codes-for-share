//#define SANITY_CHECK_REALLOCATABLE_ARRAY

#include "header.h"
#include <sys/stat.h>

int main(int argc, char* argv[]){
  PS::Initialize(argc, argv);
  
  PS::ParticleSystem<RealPtcl> sph_system;
  PS::DomainInfo dinfo;
  PS::TreeForForceShort<RESULT::Dens, EPI::Dens, EPJ::Dens>::Gather_or_Scatter  dens_tree;
  PS::TreeForForceShort<RESULT::Gradients, EPI::Gradients, EPJ::Gradients>::Gather_or_Scatter  gradients_tree;
  PS::TreeForForceShort<RESULT::Elastic, EPI::Elastic, EPJ::Elastic>::Symmetry_or_Scatter  elast_tree;
  PS::TreeForForceLong <RESULT::Grav , EPI::Grav , EPJ::Grav >::Monopole grav_tree;
  My_SPH_initialize(sph_system,dinfo,dens_tree,gradients_tree,elast_tree,grav_tree,argv[1],argc,argv);
  char filename[]="FDPS-basalt-sphere-collision-0000.txt";
  char datadir[256];
  sprintf(datadir,"data-%s",argv[2]);
  char dataname[256];
  char finalname[256];
  int step=0,i,error;
  double t_end=1.0e5;
  int tEndExponent=(int)log10(t_end+PARAM::ACC);
  double tEndMantissa=t_end/pow(10.0,tEndExponent);
  double t_interval=1000.0;
  double t_out=((int)(getTime()/t_interval)+1)*t_interval;
  int l=(int)(getTime()/t_interval);
  char midfile[]="midfile.bin";
  FileHeader header;
  time_t start,end;
  time_t sim_start,sim_now;
  if(PS::Comm::getRank()==0){
    sim_start=time(NULL);
  }
  PS::Comm::broadcast(&(sim_start),1,0);
  if(PS::Comm::getRank()==0){
    sim_now=time(NULL);
  }
  PS::Comm::broadcast(&(sim_now),1,0);

  if(getTime()<PARAM::ACC){
    if(PS::Comm::getRank()==0){
      if(mkdir(datadir,0777)!=0){
	error=1;
	printf("cannot make directory %s ! \n",datadir);
      }
      else error=0;
    }
    PS::Comm::broadcast(&error,1,0);
    if(error==1){
      PS::Finalize(); 
      exit(1);
    }
  }

  if(getTime()<PARAM::ACC){
    sprintf(filename,"FDPS-basalt-sphere-collision-%04d.txt",l);
    sprintf(dataname,"%s/%s",datadir,filename);
    header.time = getTime();
    header.Nbody = sph_system.getNumberOfParticleGlobal();
    sph_system.writeParticleAscii(dataname,header);
  }
  l++;
  
  for( ; getTime()<t_end&&(sim_now-sim_start)<60*60*11.75 ; step++){
    start=time(NULL);
    one_time_development_roop(sph_system,dinfo,dens_tree,gradients_tree,elast_tree,grav_tree);
    end=time(NULL);
    
    if(PS::Comm::getRank()==0){
      printf("roop %d is done! t=%e dt=%e %ld sec \n",step,getTime(),getTimeStep(),end-start);
    }
    
    if(getTime()>t_out){
      sprintf(filename,"FDPS-basalt-sphere-collision-%04d.txt",l);
      sprintf(dataname,"%s/%s",datadir,filename);
      header.time = getTime();
      header.Nbody = sph_system.getNumberOfParticleGlobal();
      sph_system.writeParticleAscii(dataname,header);
      l++;
      t_out += t_interval;
    }

    if(PS::Comm::getRank()==0){
      sim_now=time(NULL);
    }
    PS::Comm::broadcast(&(sim_now),1,0);

  }

  header.time = getTime();
  header.Nbody = sph_system.getNumberOfParticleGlobal();
  sph_system.writeParticleBinary(midfile,header);

  if(getTime()>=t_end){
    sprintf(finalname,"%s/final-%s-t=%2.1fe%d.bin",datadir,argv[2],tEndMantissa,tEndExponent);
    header.time = getTime();
    header.Nbody = sph_system.getNumberOfParticleGlobal();
    sph_system.writeParticleBinary(finalname,header);
  }

  free_aneos_array();
  PS::Finalize();
  return(0);
}
