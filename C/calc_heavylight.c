#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <stdint.h>
#include <errno.h>

#define Ns 32 
#define Nt 28
#define START 4000
#define STOP 4100
#define STEP 10

// index the flat prop_l/h array
size_t idx_light(int x, int y, int z, int s1, int c1, int s2, int c2, int t) {
    return (((((((size_t)x * Ns + y) * Ns + z) * 4 + s1) * 3 + c1) * 4 + s2) * 3 + c2) *Nt+ t;
}
size_t idx_heavy(int x, int y, int z, int s1, int c1, int s2, int c2, int t) {
    return (((((((size_t)x * Ns + y) * Ns + z) * 2 + s1) * 3 + c1) * 2 + s2) * 3 + c2) *Nt+t;
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
    for (int t=0;t<Nt;t++){	
    for (int x = 0; x < Ns; x++){
    	for (int y = 0; y < Ns; y++){
    	    for (int z = 0; z < Ns; z++){
    		for (int c1 = 0; c1 < 3; c1++){
    		    for (int c2 = 0; c2 < 3; c2++){
			// Extract 4x4 block for fixed x,y,z,c1,c2
			double complex block[4][4];
			for (int s1 = 0; s1 < 4; s1++){
			    for (int s2 = 0; s2 < 4; s2++){
			    	size_t idx = idx_light(x, y, z, s1, c1, s2, c2, t);
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
			    	size_t idx = idx_light(x, y, z, s1, c1, s2, c2, t);
			    	prop_tr[idx] = temp2[s1][s2] * 0.5;
		            }
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
    for (int t=0; t<Nt; t++){
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
				    fprintf(stderr, "Error reading data at t=%d x=%d y=%d z=%d s1=%d c1=%d s2=%d c2=%d\n",
                                            t, x, y, z, s1, c1, s2, c2);
				    fclose(f);
				    exit(1);
				}
				size_t idx = idx_light(x, y, z, s1, c1, s2, c2, t); 
				prop_l[idx] = repart+impart*I;
    			    }
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
    for (int t=0; t<Nt; t++){
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
				    fprintf(stderr, "Error reading data at t=%d x=%d y=%d z=%d s1=%d c1=%d s2=%d c2=%d\n",
                                           t, x, y, z, s1, c1, s2, c2);
                                    fclose(f);
                                    exit(1);
                            	}
				size_t idx = idx_heavy(x, y, z, s1, c1, s2, c2, t);
				prop_h[idx] = repart + impart*I;
			    }
                        }
                    }
                }
            }
        }
    }
    }
    fclose(f);
}

void calc_pseudo(const double complex *arrayL, const double complex *arrayH, double *pseudo){
    for (int t = 0; t < Nt; t++){
	double pscalar = 0.0;
	for (int x = 0; x < Ns; x++){
	    for (int y = 0; y < Ns; y++){
		for (int z = 0; z < Ns; z++){
		    for (int c1 = 0; c1 < 3; c1++){
			for (int c2 = 0; c2 < 3; c2++){
			// light: [x][y][z][0..1][c1][0..1][c2]
		        // heavy: [x][y][z][0..1][c1][0..1][c2]

	         	double complex l00 = arrayL[idx_light(x, y, z, 0, c1, 0, c2, t)];
			double complex l01 = arrayL[idx_light(x, y, z, 0, c1, 1, c2, t)];
			double complex l10 = arrayL[idx_light(x, y, z, 1, c1, 0, c2, t)];
			double complex l11 = arrayL[idx_light(x, y, z, 1, c1, 1, c2, t)];

			double complex h00 = arrayH[idx_heavy(x, y, z, 0, c1, 0, c2, t)];
			double complex h01 = arrayH[idx_heavy(x, y, z, 0, c1, 1, c2, t)];
			double complex h10 = arrayH[idx_heavy(x, y, z, 1, c1, 0, c2, t)];
			double complex h11 = arrayH[idx_heavy(x, y, z, 1, c1, 1, c2, t)];

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
	pseudo[t] = pscalar;
    }
}

void calc_vector(const double complex *arrayL, const double complex *arrayH, double *vector){
    
    double complex S15[2][2]={
	    {0.0 + 0.0*I, 0.0 + 1.0*I},
	    {0.0 + 1.0*I, 0.0 + 0.0*I}
    };
    double complex S25[2][2]={
            {0.0 + 0.0*I, 1.0 + 0.0*I},
            {-1.0 + 0.0*I, 0.0 + 0.0*I}
    };
    double complex S35[2][2]={
            {0.0 + 1.0*I, 0.0 + 0.0*I},
            {0.0 + 0.0*I, 0.0 - 1.0*I}
    };
    for (int t = 0; t < Nt; t++){
        double vectorval = 0.0;
        for (int x = 0; x < Ns; x++){
            for (int y = 0; y < Ns; y++){
                for (int z = 0; z < Ns; z++){
                    for (int c1 = 0; c1 < 3; c1++){
                        for (int c2 = 0; c2 < 3; c2++){
                        // light: [x][y][z][0..1][c1][0..1][c2]
                        // heavy: [x][y][z][0..1][c1][0..1][c2]

                        double complex l00 = arrayL[idx_light(x, y, z, 0, c1, 0, c2, t)];
                        double complex l01 = arrayL[idx_light(x, y, z, 0, c1, 1, c2, t)];
                        double complex l10 = arrayL[idx_light(x, y, z, 1, c1, 0, c2, t)];
                        double complex l11 = arrayL[idx_light(x, y, z, 1, c1, 1, c2, t)];

                        double complex h00 = arrayH[idx_heavy(x, y, z, 0, c1, 0, c2, t)];
                        double complex h01 = arrayH[idx_heavy(x, y, z, 0, c1, 1, c2, t)];
                        double complex h10 = arrayH[idx_heavy(x, y, z, 1, c1, 0, c2, t)];
                        double complex h11 = arrayH[idx_heavy(x, y, z, 1, c1, 1, c2, t)];

                        // Matrix product and Hermitian conjugate (transpose + conjugate)
                        // Ti = Si5.H * h * Si5 l.H, where l and h are 2x2 matrices
                        // T1/2[0,0] =  
                        double complex T100 = -S15[0][1]*h10*S15[0][1]*conj(l01) -S15[0][1]*h11*S15[1][0]*l00;
			double complex T111 = -S15[1][0]*h00*S15[0][1]*l11 -S15[1][0]*h00*S15[1][0]*conj(l10);
                        vectorval += creal(T100 + T111); //trace of T1
                        
			double complex T200 = -S25[0][1]*h10*S25[0][1]*conj(l01) -S25[0][1]*h11*S25[1][0]*l00;
                        double complex T211 = -S25[1][0]*h00*S25[0][1]*l11 -S25[1][0]*h00*S25[1][0]*conj(l10);
                        vectorval += creal(T200 + T211); //trace of T1
			
			double complex T300 = -S35[0][0]*h00*S35[0][0]*l00 -S35[0][0]*h01*S35[1][1]*conj(l01);
                        double complex T311 = -S35[1][1]*h10*S35[0][0]*conj(l10) -S35[1][1]*h11*S35[1][1]*l11;
                        vectorval += creal(T300 + T311); //trace of T1

			}
                    }
                }
            }
        }
        vector[t] = vectorval;
    }

}
//////////////////////////////////////////////////////////////////////
int main(){

double pseudo_sum[Nt];
double pseudo_avg[Nt];
double vector_sum[Nt];

//double pseudo[Nt];
//double vector[Nt];
//for (int i=0; i<Nt; i++){
//    pseudo[i] = 0.0;
//    vector[i] = 0.0;
//}
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
//for (int nt=0; nt<Nt; nt++){
    double complex *prop_l = malloc(sizeof(double complex)*Ns*Ns*Ns*4*3*4*3*Nt);
    double complex *prop_h = malloc(sizeof(double complex)*Ns*Ns*Ns*2*3*2*3*Nt);
    double complex *prop_tr = malloc(sizeof(double complex)*Ns*Ns*Ns*4*3*4*3*Nt);
    if (!prop_l || !prop_h || !prop_tr) {
    	perror("malloc prop_l/h/tr");
	free(prop_l);
	free(prop_h);
	free(prop_tr);
    	exit(1);
    }
    
    char hfile[256];
    char lfile[256];
    int ncfg = 0;
    for (int cfg=START; cfg<=STOP; cfg+=STEP){
    	double pseudo[Nt];
    	double vector[Nt];
    	for (int i=0; i<Nt; i++){
            pseudo[i] = 0.0;
            vector[i] = 0.0;
    	}
 
    	sprintf(hfile,"/data/rachel/Bmeson/nrqcd/propagators/%dx%d/prop/sprop/sprop.%dx%d_%d", 
		    Ns,Nt,Ns,Nt,cfg);
    	sprintf(lfile, "/data/rachel/Bmeson/rqcd/propagators/%dx%d/light/s0/Gen2_%dx%dcfgn%d.s0.m0", 
		    Ns,Nt,Ns,Nt,cfg);
	
	// Check heavy file
        FILE *fh = fopen(hfile, "rb");
        if (!fh) {
            printf("Warning: Cannot open heavy file for config %d, skipping.\n", cfg);
            continue;
        }
        fclose(fh);

        // Check light file
        FILE *fl = fopen(lfile, "rb");
        if (!fl) {
            printf("Warning: Cannot open light file for config %d, skipping.\n", cfg);
            continue;
        }
        fclose(fl);

    	read_heavy_prop(hfile, prop_h);
    	read_light_prop(lfile, prop_l);
    	printf("Files read successfully!\n");
    
    	transform_light(prop_l, prop_tr, trans, trans_inv);
    	calc_pseudo(prop_tr, prop_h, pseudo);
 	printf("Calculated pseudoscalar:");
    	for (int t = 0; t < Nt; t++) {
	    printf("pseudo[%d] = %f\n", t, pseudo[t]);
    	    pseudo_sum[t] += pseudo[t];
	}
	calc_vector(prop_tr, prop_h, vector);
	printf("Calculated vector:");
	for (int t = 0; t < Nt; t++) {
            printf("vector[%d] = %f\n", t, vector[t]);
            vector_sum[t] += vector[t];
        }
	ncfg++;
    }
    
    for (int t=0;t <=Nt; t++){
	pseudo_avg[t] = pseudo_sum[t] / ncfg;
    }
    printf("\nAveraged pseudoscalar over %d configs:\n", ncfg);
    for (int t = 0; t < Nt; t++) {
        printf("pseudo_avg[%d] = %f\n", t, pseudo_avg[t]);
    }

    free(prop_l);
    free(prop_h);
    free(prop_tr);


    return 0;

}

