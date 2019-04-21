#ifndef _IO_JING_H_
#define _IO_JING_H_
#include <stdint.h>
#include "../particle.h"
 
struct jing_header{
	uint32_t np;
	uint32_t ips;
	float ztp;
	float omegat;
	float lamdat;
	float boxsize;
	float xscale;
	float vscale;
};

struct ljing_header{
	int64_t np;
	int64_t ips;
	float ztp;
	float omegat;
	float lamdat;
	float boxsize;
	float xscale;
	float vscale;
};

void load_particles_jing(char *filename, struct particle **p, int64_t *num_p);
void load_particles_ljing(char *filename, struct particle **p, int64_t *num_p);

#endif /* _IO_JING_H_ */
