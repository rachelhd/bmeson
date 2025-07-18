#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <stdint.h>
#include <errno.h>

#define Ns 32 
#define Nt 28

// index the flat prop_l/h array
size_t idx_light(int x, int y, int z, int s1, int c1, int s2, int c2) {
    return ((((((size_t)x * Ns + y) * Ns + z) * 4 + s1) * 3 + c1) * 4 + s2) * 3 + c2;
}
size_t idx_heavy(int x, int y, int z, int s1, int c1, int s2, int c2) {
    return ((((((size_t)x * Ns + y) * Ns + z) * 2 + s1) * 3 + c1) * 2 + s2) * 3 + c2;
}

//read binary light prop file
void read_light_prop(const char *filename, 
                     double complex *prop_l){
    FILE *f = fopen(filename, "rb");
    if (!f){
        //fprintf(stderr, "Unable to open file '%s': %d\n", filename, strerror(errno));
        perror("Error opening file");
	exit(1);
    }
    
    // Read header
    int header[4];
    double norm;
    if (fread(header,sizeof(int),4,f)!=4 ||
    	fread(&norm,sizeof(double),1,f)!=1) {
    	fclose(f);
	exit(1);	
    }

    //Read prop data
    for (int x=0; x<Ns; x++){
  	for (int y=0; y<Ns; y++){
    	    for (int z=0; z<Ns; z++){
    		for (int s1=0; s1<4; s1++){
		    for (int c1=0; c1<3; c1++){
    			for (int s2=0; s2<4; s2++){
			    for (int c2=0; c2<3; c2++){
    				double repart, impart;
				if (fread(&repart,sizeof(double),1,f)!=1 ||
				    fread(&impart,sizeof(double),1,f)!=1){
				    fprintf(stderr, "Error reading data at x=%d y=%d z=%d s1=%d c1=%d s2=%d c2=%d\n",
                                            x, y, z, s1, c1, s2, c2);
				    fclose(f);
				    exit(1);
				}
				size_t idx = idx_light(x, y, z, s1, c1, s2, c2); 
				prop_l[idx] = repart+impart*I;
    			    }
    			}
		    }
		}
	    }
	}
    }
    fclose(f);
}

//read binary heavy prop file
void read_heavy_prop(const char *filename,
                     double complex *prop_h){
    FILE *f = fopen(filename, "rb");
    if (!f){
        //fprintf(stderr, "Unable to open file '%s': %d\n", filename, strerror(errno));
        perror("Error opening file");
	exit(1);
        }

    // Read header
    int header;
    if (fread(&header,sizeof(int),1,f) !=1){
	perror("Error reading header from heavy file");	
    } 

    //Read prop data
    for (int s1=0; s1<2; s1++){
        for (int s2=0; s2<2; s2++){
            for (int c1=0; c1<3; c1++){
                for (int c2=0; c2<3; c2++){
                    for (int z=0; z<Ns; z++){
                        for (int y=0; y<Ns; y++){
                            for (int x=0; x<Ns; x++){
                                double repart, impart;
                                if (fread(&repart,sizeof(double),1,f)!=1 ||
                                    fread(&impart,sizeof(double),1,f)!=1){
				    fprintf(stderr, "Error reading data at x=%d y=%d z=%d s1=%d c1=%d s2=%d c2=%d\n",
                                            x, y, z, s1, c1, s2, c2);
                                    fclose(f);
                                    exit(1);
                            	}
				size_t idx = idx_heavy(x, y, z, s1, c1, s2, c2);
				prop_h[idx] = repart + impart*I;
			    }
                        }
                    }
                }
            }
        }
    }
    fclose(f);
}


int main(){

//double complex prop_l[Ns][Ns][Ns][4][3][4][3];
//double complex prop_h[Ns][Ns][Ns][2][3][2][3];
double complex *prop_l = malloc(sizeof(double complex)*Ns*Ns*Ns*4*3*4*3);
double complex *prop_h = malloc(sizeof(double complex)*Ns*Ns*Ns*2*3*2*3);
if (!prop_l) {
    perror("malloc prop_l");
    exit(1);
}
if (!prop_h) {
    perror("malloc prop_h");
    exit(1);
}
double pseudo[Nt];
double vector[Nt];

char hfile[256];
char lfile[256];
int i = 3510;
sprintf(hfile, "/data/rachel/Bmeson/nrqcd/propagators/%dx%d/prop/sprop/sprop.%dx%d_%d", Ns,Nt,Ns,Nt,i);
sprintf(lfile, "/data/rachel/Bmeson/rqcd/propagators/%dx%d/light/s0/Gen2_%dx%dcfgn%d.s0.m0", Ns,Nt,Ns,Nt,i);

read_heavy_prop(hfile, prop_h);
read_light_prop(lfile, prop_l);
printf("Files read successfully!\n");

free(prop_l);
free(prop_h);
return 0;

}

