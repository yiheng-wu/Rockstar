#ifdef ENABLE_HDF5
#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h> /* HDF5 required */
#include <inttypes.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include "io_hdf5.h"
#include "io_tng.h"
#include "io_util.h"
#include "../universal_constants.h"
#include "../check_syscalls.h"
#include "../config_vars.h"
#include "../config.h"
#include "../particle.h"

void load_paticles_hdf5(char *filename, struct particle **p,int64_t *num_p)
{
	hid_t fid=check_H5Fopen(filename,H5F_ACC_RDONLY);

	char *dm="PartType1";
	char *loc="PartType1/Coordinates";
	char *vel="PartType1/Velocities";
	char *id="PartType1/ParticleIDs";
	herr_t status = H5Lget_info(fid,dm, NULL, H5P_DEFAULT);
	if (status == 0)
	{
		int64_t i,j;
		struct position{ double pos[3];}; 
		struct position *ploc;
		struct velocity{ float vel[3];}; 
		struct velocity *pvel;
		struct identity{ int64_t id;}; 
		struct identity *pid;

		hid_t locid=check_H5Dopen(fid,loc,dm,filename);
		hid_t velid=check_H5Dopen(fid,vel,dm,filename);
		hid_t idid=check_H5Dopen(fid,id,dm,filename);

		hid_t locspace=check_H5Dget_space(locid);
		//hid_t velspace=check_H5Dget_space(velid);
		//hid_t idspace=check_H5Dget_space(idid);

		hsize_t *dm_particle;
		status=H5Sget_simple_extent_dims(locspace,dm_particle,NULL);		
  		*p = (struct particle *)check_realloc(*p, ((*num_p)+(*dm_particle))*sizeof(struct particle), "Allocating particles.");

		ploc=(struct position *)malloc((*dm_particle)*sizeof(struct position));
		check_H5Dread(locid, H5T_NATIVE_DOUBLE,ploc,loc,dm,filename);
		H5Dclose(locid);
		for(i=0;i<(*dm_particle);i++)
		{
			for(j=0;j<3;j++)
			{
				(*p)[*num_p+i].pos[j]=ploc[i].pos[j];
			}
		}
		free(ploc);
		
		pvel=(struct velocity *)malloc((*dm_particle)*sizeof(struct velocity));
		check_H5Dread(velid, H5T_NATIVE_DOUBLE,pvel,vel,dm,filename);
		H5Dclose(velid);
		for(i=0;i<(*dm_particle);i++)
		{
			for(j=0;j<3;j++)
			{
				(*p)[*num_p+i].pos[j+3]=pvel[i].vel[j];
			}
		}
		free(pvel);

		pid=(struct identity *)malloc((*dm_particle)*sizeof(struct identity));
		check_H5Dread(idid, H5T_NATIVE_DOUBLE,pid,id,dm,filename);
		H5Dclose(idid);
		for(i=0;i<(*dm_particle);i++)
		{
			(*p)[*num_p+i].id=pid[i].id;
		}
		free(pid);

		*num_p+=*dm_particle;
	}
	else
		printf("The group either does NOT exist or some other error occurred.\n");
	H5Fclose(fid);
}

#endif /* ENABLE_HDF5 */
