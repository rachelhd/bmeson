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

void mat4_mult(const double complex A[4][4], const double complex B[4][4], double complex C[4][4]) {
    for (int i = 0; i < 4; i++)
    	for (int j = 0; j < 4; j++) {
            C[i][j] = 0.0;
            for (int k = 0; k < 4; k++){
          	C[i][j] += A[i][k] * B[k][j];
      	    }
	}
    }

// transform light to NR rep
void transform_light(
    const double complex *prop_l,  // input, chiral basis
    double complex *prop_tr,     // output, NR basis
    const double complex trans[4][4], //transform mat
    const double complex trans_inv[4][4] 
) {
    for (int x = 0; x < Ns; x++){
    	for (int y = 0; y < Ns; y++){
    	    for (int z = 0; z < Ns; z++){
    		for (int c1 = 0; c1 < 3; c1++){
    		    for (int c2 = 0; c2 < 3; c2++){
			// Extract 4x4 block for fixed x,y,z,c1,c2
			double complex block[4][4];
			for (int s1 = 0; s1 < 4; s1++){
			    for (int s2 = 0; s2 < 4; s2++){
			    	size_t idx = idx_light(x, y, z, s1, c1, s2, c2);
			    	block[s1][s2] = prop_l[idx];
			    }
			}
			// temp1 = trans * block
			double complex temp1[4][4];
			mat4_mult(trans, block, temp1);

			// temp2 = temp1 * trans_inv
			double complex temp2[4][4];
			mat4_mult(temp1, trans_inv, temp2);

			// Store result back, scaled by 1/2
			for (int s1 = 0; s1 < 4; s1++){
			    for (int s2 = 0; s2 < 4; s2++){
			    	size_t idx = idx_light(x, y, z, s1, c1, s2, c2);
			    	prop_tr[idx] = temp2[s1][s2] * 0.5;
		            }
			}
		    }
		}
	    }
	}
    }
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
	fprintf(stderr, "Error reading header/norm from file %s\n", filename);
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
	fclose(f);
	exit(1);	
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

double calc_pseudo(const double complex *arrayL, const double complex *arrayH) {
    double pscalar = 0.0;
    for (int x = 0; x < Ns; x++){
    	for (int y = 0; y < Ns; y++){
    	    for (int z = 0; z < Ns; z++){
    		for (int c1 = 0; c1 < 3; c1++){
    		    for (int c2 = 0; c2 < 3; c2++){
        	    // light: [x][y][z][0..1][c1][0..1][c2]
        	    // heavy: [x][y][z][0..1][c1][0..1][c2]

		    double complex l00 = arrayL[idx_light(x, y, z, 0, c1, 0, c2)];
	 	    double complex l01 = arrayL[idx_light(x, y, z, 0, c1, 1, c2)];
		    double complex l10 = arrayL[idx_light(x, y, z, 1, c1, 0, c2)];
		    double complex l11 = arrayL[idx_light(x, y, z, 1, c1, 1, c2)];

		    double complex h00 = arrayH[idx_heavy(x, y, z, 0, c1, 0, c2)];
		    double complex h01 = arrayH[idx_heavy(x, y, z, 0, c1, 1, c2)];
		    double complex h10 = arrayH[idx_heavy(x, y, z, 1, c1, 0, c2)];
		    double complex h11 = arrayH[idx_heavy(x, y, z, 1, c1, 1, c2)];

		    // Matrix product and Hermitian conjugate (transpose + conjugate)
		    // S = l * h.H, where l and h are 2x2 matrices
		    // S[0,0] = l00*conj(h00) + l01*conj(h01)
		    // S[1,1] = l10*conj(h10) + l11*conj(h11)
		    double complex S00 = l00*conj(h00) + l01*conj(h01);
		    double complex S11 = l10*conj(h10) + l11*conj(h11);

		    pscalar += creal(S00 + S11); //trace of S
	            }
    		}
	    }
	}
    }
    return pscalar;
}

//////////////////////////////////////////////////////////////////////
int main(){

double pseudo[Nt];
double vector[Nt];
for (int i=0; i<Nt; i++){
    pseudo[i] = 0.0;
    vector[i] = 0.0;
}
double complex trans[4][4] = {
    {1.0 + 0.0*I, 0.0 + 0.0*I, -1.0 + 0.0*I, 0.0 + 0.0*I},
    {0.0 + 0.0*I, 1.0 + 0.0*I,  0.0 + 0.0*I, -1.0 + 0.0*I},
    {1.0 + 0.0*I, 0.0 + 0.0*I,  1.0 + 0.0*I,  0.0 + 0.0*I},
    {0.0 + 0.0*I, 1.0 + 0.0*I,  0.0 + 0.0*I,  1.0 + 0.0*I}
};
double complex trans_inv[4][4];
for (int i = 0; i < 4; i++){
    for (int j = 0; j < 4; j++){
        trans_inv[i][j] = conj(trans[j][i]);
    }
}
for (int nt=0; nt<Nt; nt++){
    double complex *prop_l = malloc(sizeof(double complex)*Ns*Ns*Ns*4*3*4*3);
    double complex *prop_h = malloc(sizeof(double complex)*Ns*Ns*Ns*2*3*2*3);
    double complex *prop_tr = malloc(sizeof(double complex)*Ns*Ns*Ns*4*3*4*3);
    if (!prop_l || !prop_h || !prop_tr) {
    	perror("malloc prop_l/h/tr");
	free(prop_l);
	free(prop_h);
	free(prop_tr);
    	exit(1);
    }
    
    char hfile[256];
    char lfile[256];
    int i = 3510;
    sprintf(hfile, "/data/rachel/Bmeson/nrqcd/propagators/%dx%d/prop/sprop/sprop.%dx%d_%d", 
		    Ns,Nt,Ns,Nt,i);
    sprintf(lfile, "/data/rachel/Bmeson/rqcd/propagators/%dx%d/light/s0/Gen2_%dx%dcfgn%d.s0.m0", 
		    Ns,Nt,Ns,Nt,i);

    read_heavy_prop(hfile, prop_h);
    read_light_prop(lfile, prop_l);
    if (nt == Nt-1){
        printf("Files read successfully!\n");
    }
    transform_light(prop_l, prop_tr, trans, trans_inv);
    pseudo[nt] = calc_pseudo(prop_tr, prop_h);
    if (nt == Nt-1){
    	printf("Calculated pseudoscalar:");
	for (int t = 0; t < Nt; t++) {
	    printf("pseudo[%d] = %f\n,", t, pseudo[t]);
	}
    }
    free(prop_l);
    free(prop_h);
    free(prop_tr);
}

return 0;

}

