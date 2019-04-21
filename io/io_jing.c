#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include "io_jing.h"
#include "io_util.h"
#include "unistd.h"
#include "../universal_constants.h"
#include "../check_syscalls.h"
#include "../config_vars.h"
#include "../config.h"
#include "../particle.h"


int64_t NF=2147483639;

//Parameter from Jing simulation 6411
float alpha=1,nstepf=600,dp=0.12;


int dummy;
#define SKIP(fp) fread(&dummy,sizeof(dummy),1,fp)
//int大小端转换
unsigned int UISWAP(int x)
{
	return \
	(x&0xff000000)>>24|\
	(x&0x00ff0000)>>8|\
	(x&0x0000ff00)<<8|\
	(x&0x000000ff)<<24;
}
//float大小端转换
typedef union FLOAT_CHAR
{
	float f;
	char c[4];
}float_char;
float FSWAP(float x)
{
	float_char f1,f2;
	f1.f=x;
	f2.c[0]=f1.c[3];
	f2.c[1]=f1.c[2];
	f2.c[2]=f1.c[1];
	f2.c[3]=f1.c[0];
	return f2.f;
}
struct jing_header jing_extract_header_info(FILE *pos){
	struct jing_header header; 
	//read header infomation from position file
	SKIP(pos);
	fread(&header,sizeof(header),1,pos);
	SKIP(pos);
	fprintf(stderr,"SKIP=%d\n",UISWAP(dummy));
	
	header.np=UISWAP(header.np);
	header.ips=UISWAP(header.ips);
	header.ztp=FSWAP(header.ztp);
	header.omegat=FSWAP(header.omegat);
	header.lamdat=FSWAP(header.lamdat);
	header.boxsize=FSWAP(header.boxsize);
	header.xscale=FSWAP(header.xscale);
	header.vscale=FSWAP(header.vscale);


	float scalep,scalepf,a,h,vfact2;

	scalep=header.ips*JING_DP+1;
	scalepf=JING_NSTEP*JING_DP+1;
	a=scalep/scalepf;
	h=sqrt(Om/pow(a,3)+1-Om);
	vfact2=JING_ALPHA*scalep;
	
	
	BOX_SIZE=header.boxsize;
	TOTAL_PARTICLES=header.np;
	SCALE_NOW=a;
	PARTICLE_MASS=CRITICAL_DENSITY*pow(BOX_SIZE,3)*Om/TOTAL_PARTICLES;
	AVG_PARTICLE_SPACING=cbrt(PARTICLE_MASS/(Om*CRITICAL_DENSITY));

	//to calculate velocity
	header.vscale=vfact2*BOX_SIZE*100*h*SCALE_NOW;

	return header;
}
void jing_read_id(FILE *id,struct particle *p,int64_t num_p,int64_t np){
	int64_t j,nb=0;//nb,number of bytes read already
	int32_t tem;
	SKIP(id);
	for(j=0;j<np;j++)
	{
		fread(&tem,sizeof(tem),1,id);
		p[num_p+j].id=UISWAP(tem);
		nb+=4;
		if(nb==NF-3)
		{
			float_char joint;
			fread(&joint.c[0],1,1,id);
			fread(&joint.c[1],1,1,id);
			fread(&joint.c[2],1,1,id);
			SKIP(id);
			SKIP(id);
			fread(&joint.c[3],1,1,id);
			p[num_p+j+1].id=UISWAP(joint.f);
			j++;
			nb+=4;
		}
		else if(nb==NF*2-2)
		{
			float_char joint;
			fread(&joint.c[0],1,1,id);
			fread(&joint.c[1],1,1,id);
			SKIP(id);
			SKIP(id);
			fread(&joint.c[2],1,1,id);
			fread(&joint.c[3],1,1,id);
			p[num_p+j+1].id=UISWAP(joint.f);
			j++;
			nb+=4;
		}
		else if(nb==NF*3-1)
		{
			float_char joint;
			fread(&joint.c[0],1,1,id);
			SKIP(id);
			SKIP(id);
			fread(&joint.c[1],1,1,id);
			fread(&joint.c[2],1,1,id);
			fread(&joint.c[3],1,1,id);
			p[num_p+j+1].id=UISWAP(joint.f);
			j++;
			nb+=4;
		}
		else if(nb==NF*4)
		{
			SKIP(id);
			SKIP(id);
			nb=0;
		}
	}
}
void jing_read_pos(FILE *pos,struct particle *p,int64_t num_p,struct jing_header header){
	int64_t i,j,nb=0;//nb,number of bytes read already
	int64_t np;
	float boxsize;
	np=header.np;
	boxsize=header.boxsize;

	float tem;
	SKIP(pos);
	for(i=0;i<3;i++)
	{
		for(j=0;j<np;j++)
		{
			fread(&tem,sizeof(tem),1,pos);
			p[num_p+j].pos[i]=FSWAP(tem)*boxsize;
			nb+=4;
			if(nb==NF-3)
			{
				float_char joint;
				fread(&joint.c[0],1,1,pos);
				fread(&joint.c[1],1,1,pos);
				fread(&joint.c[2],1,1,pos);
				SKIP(pos);
				SKIP(pos);
				fread(&joint.c[3],1,1,pos);
				p[num_p+j+1].pos[i]=FSWAP(joint.f)*boxsize;
				j++;
				nb+=4;
			}
			else if(nb==NF*2-2)
			{
				float_char joint;
				fread(&joint.c[0],1,1,pos);
				fread(&joint.c[1],1,1,pos);
				SKIP(pos);
				SKIP(pos);
				fread(&joint.c[2],1,1,pos);
				fread(&joint.c[3],1,1,pos);
				p[num_p+j+1].pos[i]=FSWAP(joint.f)*boxsize;
				j++;
				nb+=4;
			}
			else if(nb==NF*3-1)
			{
				float_char joint;
				fread(&joint.c[0],1,1,pos);
				SKIP(pos);
				SKIP(pos);
				fread(&joint.c[1],1,1,pos);
				fread(&joint.c[2],1,1,pos);
				fread(&joint.c[3],1,1,pos);
				p[num_p+j+1].pos[i]=FSWAP(joint.f)*boxsize;
				j++;
				nb+=4;
			}
			else if(nb==NF*4)
			{
				SKIP(pos);
				SKIP(pos);
				nb=0;
			}
		}
	}
}
void jing_read_vel(FILE *vel,struct particle *p,int64_t num_p,struct jing_header header){
	int64_t i,j,nb=0;//nb,number of bytes read already
	float tem;
	int64_t np;
	float vscale;
	np=header.np;
	vscale=header.vscale;
	
	struct jing_header header_dummy;
	SKIP(vel);
	fread(&header_dummy,sizeof(header_dummy),1,vel);
	SKIP(vel);
	SKIP(vel);
	for(i=3;i<6;i++)
	{
		for(j=0;j<np;j++)
		{
			fread(&tem,sizeof(tem),1,vel);
			p[num_p+j].pos[i]=FSWAP(tem)*vscale;
//	fprintf(stderr,"xyz=%f\n",p[num_p+j].pos[i]);
			nb+=4;
			if(nb==NF-3)
			{
				float_char joint;
				fread(&joint.c[0],1,1,vel);
				fread(&joint.c[1],1,1,vel);
				fread(&joint.c[2],1,1,vel);
				SKIP(vel);
				SKIP(vel);
				fread(&joint.c[3],1,1,vel);
				p[num_p+j+1].pos[i]=FSWAP(joint.f)*vscale;
				j++;
				nb+=4;
			}
			else if(nb==NF*2-2)
			{
				float_char joint;
				fread(&joint.c[0],1,1,vel);
				fread(&joint.c[1],1,1,vel);
				SKIP(vel);
				SKIP(vel);
				fread(&joint.c[2],1,1,vel);
				fread(&joint.c[3],1,1,vel);
				p[num_p+j+1].pos[i]=FSWAP(joint.f)*vscale;
				j++;
				nb+=4;
			}
			else if(nb==NF*3-1)
			{
				float_char joint;
				fread(&joint.c[0],1,1,vel);
				SKIP(vel);
				SKIP(vel);
				fread(&joint.c[1],1,1,vel);
				fread(&joint.c[2],1,1,vel);
				fread(&joint.c[3],1,1,vel);
				p[num_p+j+1].pos[i]=FSWAP(joint.f)*vscale;
				j++;
				nb+=4;
			}
			else if(nb==NF*4)
			{
				SKIP(vel);
				SKIP(vel);
				nb=0;
			}
		}
	}
}

void load_particles_jing(char *filename, struct particle **p, int64_t *num_p)
{

	struct jing_header header;
	char id_fn[1024],pos_fn[1024],vel_fn[1024];//fn stands for filename

	int i,j;
	char * compare="/";
	for(i=0;i<1023;i++)
	{
		if (filename[i]==compare[0] && filename[i+1]==compare[0])
		{
			for(j=0;j<i+2;j++)
			{
				id_fn[j]=filename[j];
				pos_fn[j]=filename[j];
				vel_fn[j]=filename[j];
			}
			sprintf(id_fn,"%sid",id_fn);
			sprintf(pos_fn,"%spos",pos_fn);
			sprintf(vel_fn,"%svel",vel_fn);
			for(j=i+2;j<1021;j++)
			{
				id_fn[j+2]=filename[j];
				pos_fn[j+3]=filename[j];
				vel_fn[j+3]=filename[j];
			}
			break;
		}
	}

	FILE *id,*pos,*vel;
	id=check_fopen(id_fn,"rb");
	pos=check_fopen(pos_fn,"rb");
	vel=check_fopen(vel_fn,"rb");


	header=jing_extract_header_info(pos);

	*p=(struct particle *)check_realloc(*p,((*num_p)+header.np)*sizeof(struct particle),"Allocating particles.");

	jing_read_vel(vel,*p,*num_p,header);
	jing_read_pos(pos,*p,*num_p,header);
	jing_read_id(id,*p,*num_p,header.np);

//	profile(*p,header.np,*num_p);

	fclose(id);
	fclose(pos);
	fclose(vel);

	*num_p+=header.np;
//	fprintf(stderr,"num_p=%lld\n",*num_p);
}

/***********************************LLLLJING**********************************/




struct ljing_header ljing_extract_header_info(FILE *pos,char *head_info){
	struct ljing_header header; 
	//read header infomation from position file
	fread(&header,sizeof(header),1,pos);
	SKIP(pos);
	SKIP(pos);
	fprintf(stderr,"SKIP=%d\n",dummy);

	FILE *h_info=fopen(head_info,"w");

	float scalep,scalepf,a,h,vfact2;

	scalep=header.ips*JING_DP+1;
	scalepf=JING_NSTEP*JING_DP+1;
	a=scalep/scalepf;
	h=sqrt(Om/pow(a,3)+1-Om);
	vfact2=JING_ALPHA*scalep;
	
	
	BOX_SIZE=header.boxsize;
	TOTAL_PARTICLES=header.np;
	SCALE_NOW=a;
	PARTICLE_MASS=CRITICAL_DENSITY*pow(BOX_SIZE,3)*Om/TOTAL_PARTICLES;
	AVG_PARTICLE_SPACING=cbrt(PARTICLE_MASS/(Om*CRITICAL_DENSITY));


	fprintf(h_info,"%f\n",BOX_SIZE);
	fprintf(h_info,"%lld\n",TOTAL_PARTICLES);
	fprintf(h_info,"%f\n",SCALE_NOW);
	fprintf(h_info,"%f\n",PARTICLE_MASS);
	fprintf(h_info,"%f\n",AVG_PARTICLE_SPACING);

	//to calculate velocity
	LJING_VSCALE=vfact2*BOX_SIZE*100*h*SCALE_NOW;
	fprintf(h_info,"%f\n",LJING_VSCALE);
	
	fclose(h_info);

	return header;
}
void ljing_read_head(char * head_info){
	FILE *h_info=check_fopen(head_info,"r");
	fscanf(h_info,"%lf\n",&BOX_SIZE);
	fscanf(h_info,"%lld\n",&TOTAL_PARTICLES);
	fscanf(h_info,"%lf\n",&SCALE_NOW);
	fscanf(h_info,"%lf\n",&PARTICLE_MASS);
	fscanf(h_info,"%lf\n",&AVG_PARTICLE_SPACING);
	fscanf(h_info,"%lf\n",&LJING_VSCALE);
	fclose(h_info);
}
void ljing_read_id(FILE *id,struct particle *p,int64_t num_p){


}
void ljing_read_pos(FILE *pos,struct particle *p,int64_t num_p){
	int64_t i,j,nb=0;//nb,number of bytes read already

	float tem;
//	SKIP(pos);
	for(j=0;j<TOTAL_PARTICLES/NUM_BLOCKS;j++)
	{
		for(i=0;i<3;i++)
		{
			fread(&tem,sizeof(tem),1,pos);
			p[num_p+j].pos[i]=tem*BOX_SIZE;
			nb+=4;
			if(nb==NF-3)
			{
				float_char joint;
				fread(&joint.c[0],1,1,pos);
				fread(&joint.c[1],1,1,pos);
				fread(&joint.c[2],1,1,pos);
				SKIP(pos);
	fprintf(stderr,"SKIP=%d\n",dummy);
				SKIP(pos);
				fread(&joint.c[3],1,1,pos);
				p[num_p+j].pos[2]=joint.f*BOX_SIZE;
				i++;
				nb+=4;
			}
			else if(nb==NF*2-2)
			{
				float_char joint;
				fread(&joint.c[0],1,1,pos);
				fread(&joint.c[1],1,1,pos);
				SKIP(pos);
	fprintf(stderr,"SKIP=%d\n",dummy);
				SKIP(pos);
				fread(&joint.c[2],1,1,pos);
				fread(&joint.c[3],1,1,pos);
				p[num_p+j].pos[2]=joint.f*BOX_SIZE;
				i++;
				nb+=4;
			}
			else if(nb==NF*3-1)
			{
				float_char joint;
				fread(&joint.c[0],1,1,pos);
				SKIP(pos);
	fprintf(stderr,"SKIP=%d\n",dummy);
				SKIP(pos);
				fread(&joint.c[1],1,1,pos);
				fread(&joint.c[2],1,1,pos);
				fread(&joint.c[3],1,1,pos);
				p[num_p+j].pos[2]=joint.f*BOX_SIZE;
				i++;
				nb+=4;
			}
			else if(nb==NF*4)
			{
				SKIP(pos);
	fprintf(stderr,"SKIP=%d\n",dummy);
				SKIP(pos);
				nb=0;
			}
		}
	}
}
void ljing_read_vel(FILE *vel,struct particle *p,int64_t num_p){
	int64_t i,j,nb=0;//nb,number of bytes read already
	struct ljing_header skip_header;
	float tem;
	SKIP(vel);
	if(dummy==40)
	{
		fread(&skip_header,sizeof(skip_header),1,vel);
		SKIP(vel);
		SKIP(vel);
	}

	fprintf(stderr,"Total particles=%lld\n",TOTAL_PARTICLES);
	for(j=0;j<TOTAL_PARTICLES/NUM_BLOCKS;j++)
	{
		for(i=3;i<6;i++)
		{
			fread(&tem,sizeof(tem),1,vel);
			p[num_p+j].pos[i]=tem*LJING_VSCALE;
	//fprintf(stderr,"xyz=%f\n",p[num_p+j].pos[i]);
			nb+=4;
			if(nb==NF-3)
			{
				float_char joint;
				fread(&joint.c[0],1,1,vel);
				fread(&joint.c[1],1,1,vel);
				fread(&joint.c[2],1,1,vel);
				SKIP(vel);
	fprintf(stderr,"VELSKIP=%d\n",dummy);
				SKIP(vel);
				fread(&joint.c[3],1,1,vel);
				p[num_p+j].pos[5]=joint.f*LJING_VSCALE;
				i++;
				nb+=4;
			}
			else if(nb==NF*2-2)
			{
				float_char joint;
				fread(&joint.c[0],1,1,vel);
				fread(&joint.c[1],1,1,vel);
				SKIP(vel);
	fprintf(stderr,"SKIP=%d\n",dummy);
				SKIP(vel);
				fread(&joint.c[2],1,1,vel);
				fread(&joint.c[3],1,1,vel);
				p[num_p+j].pos[5]=joint.f*LJING_VSCALE;
				i++;
				nb+=4;
			}
			else if(nb==NF*3-1)
			{
				float_char joint;
				fread(&joint.c[0],1,1,vel);
				SKIP(vel);
	fprintf(stderr,"SKIP=%d\n",dummy);
				SKIP(vel);
				fread(&joint.c[1],1,1,vel);
				fread(&joint.c[2],1,1,vel);
				fread(&joint.c[3],1,1,vel);
				p[num_p+j].pos[5]=joint.f*LJING_VSCALE;
				i++;
				nb+=4;
			}
			else if(nb==NF*4)
			{
				SKIP(vel);
	fprintf(stderr,"SKIP=%d\n",dummy);
				SKIP(vel);
				nb=0;
			}
		}
	}
}

void profile(struct particle *p,int np,int64_t num_p)
{
	int pro_data[1200][1200];
	int i,j;
	for(i=0;i<1200;i++)
	{
		for(j=0;j<1200;j++)
		{
			pro_data[i][j]=0;
		}
	}
	for(i=0;i<np;i++)
	{
		if (i%1000==0)
		{
			fprintf(stderr,"x=%f\n",p[num_p+i].pos[0]);
			fprintf(stderr,"y=%f\n",p[num_p+i].pos[1]);
			fprintf(stderr,"z=%f\n",p[num_p+i].pos[2]);
		}
		if (p[num_p+i].pos[0]>500 && p[num_p+i].pos[0]<500.5)
		{
			pro_data[(int)(p[num_p+i].pos[1])][(int)(p[num_p+i].pos[2])]+=1;
			//fprintf(stderr,"%lld\n",p[i].id);
		}
	}
	

	FILE *new=fopen("profile.txt","a");
	long long All=0;
	for(i=0;i<1200;i++)
	{
		for(j=0;j<1200;j++)
		{
			All+=pro_data[i][j];
			fprintf(new,"%d\n",pro_data[i][j]);
		}
	}
	fclose(new);
	printf("ALL=%lld\n",All);
}
void load_particles_ljing(char *filename, struct particle **p, int64_t *num_p)
{

	struct ljing_header header;
	char id_fn[1024],pos_fn[1024],vel_fn[1024];//fn stands for filename

	//To save the head info into temporal file for different process to read
	char head_info[200],head_tem[9];
	sprintf(head_info,"%s",OUTBASE);

	int i,j;
	char * compare="/";
	for(i=0;i<1023;i++)
	{
		if (filename[i]==compare[0] && filename[i+1]==compare[0])
		{
			for(j=0;j<i+2;j++)
			{
				id_fn[j]=filename[j];
				pos_fn[j]=filename[j];
				vel_fn[j]=filename[j];
			}
			sprintf(id_fn,"%sid",id_fn);
			sprintf(pos_fn,"%spos",pos_fn);
			sprintf(vel_fn,"%svel",vel_fn);
			for(j=i+2;j<1021;j++)
			{
				id_fn[j+2]=filename[j];
				pos_fn[j+3]=filename[j];
				vel_fn[j+3]=filename[j];
			}
			for(j=0;j<9;j++)
			{
				head_tem[j]=filename[i+j+2];
			}
			break;
		}
	}
	sprintf(head_info,"%s%s",head_info,head_tem);

	FILE *id,*pos,*vel;

	if(!IGNORE_PARTICLE_IDS) id=check_fopen(id_fn,"rb");
	pos=check_fopen(pos_fn,"rb");
	vel=check_fopen(vel_fn,"rb");

	SKIP(pos);//判断是否有header
	if(dummy==40) header=ljing_extract_header_info(pos,head_info);
	else
	{
		sleep(10);//等待含header信息的block被读取，并写入文件供其他进程读取
		ljing_read_head(head_info);
	}

	*p=(struct particle *)check_realloc(*p,((*num_p)+TOTAL_PARTICLES/NUM_BLOCKS)*sizeof(struct particle),"Allocating particles.");

	if(!IGNORE_PARTICLE_IDS) ljing_read_id(id,*p,*num_p);
	ljing_read_pos(pos,*p,*num_p);
	profile(*p,TOTAL_PARTICLES/NUM_BLOCKS,*num_p);
	ljing_read_vel(vel,*p,*num_p);


	if(!IGNORE_PARTICLE_IDS) fclose(id);
	fclose(pos);
	fclose(vel);
	*num_p+=TOTAL_PARTICLES/NUM_BLOCKS;
}























































/*
  FILE *input;
  char tag[10] = {0};
  struct gadget_header header;
  int64_t total_particles, skip, skip2, halo_particles, i;

  input = check_fopen(filename, "rb");
  gadget2_detect_filetype(input, filename);

#define gadget_variant_block(a)		\
  if (GADGET_VARIANT) { \
    fread_fortran(tag, sizeof(char)*4*GADGET_VARIANT,1,input, SWAP_ENDIANNESS);\
    if (SWAP_ENDIANNESS) swap_endian_4byte((int8_t *)(&tag)); \
    assert(!strncmp(tag, a, strlen(a)));		      \
  }

  gadget_variant_block("HEAD");
  fread_fortran(&header, GADGET_HEADER_SIZE, 1, input, SWAP_ENDIANNESS);
  gadget2_extract_header_info(&header);

  total_particles = 0;
  for (i=0; i<6; i++) total_particles += header.num_particles[i];

  if (GADGET_SKIP_NON_HALO_PARTICLES) {
    for (skip=0,i=0; i<GHPT; i++) skip += header.num_particles[i];
    halo_particles = header.num_particles[GHPT];
    skip2 = total_particles - skip - halo_particles;
  } else {
    skip = skip2 = 0;
    halo_particles = total_particles;
  }


  gadget_variant_block("POS");
  gadget2_read_stride(input, *num_p, halo_particles, 3, sizeof(float), 
    *p, (char *)&(p[0][0].pos[0])-(char*)(p[0]), skip, skip2, filename);

  gadget_variant_block("VEL");
  gadget2_read_stride(input, *num_p, halo_particles, 3, sizeof(float), 
    *p, (char *)&(p[0][0].pos[3])-(char*)(p[0]), skip, skip2, filename);

  gadget_variant_block("ID");
  gadget2_read_stride(input, *num_p, halo_particles, 1, GADGET_ID_BYTES, 
    *p, (char *)&(p[0][0].id)-(char*)(p[0]), skip, skip2, filename);

  gadget2_rescale_particles(*p, *num_p, halo_particles);

  *num_p += halo_particles;
  fclose(input);
}

*/


