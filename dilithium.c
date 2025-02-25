#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<gmp.h>
#include<math.h>
#include<stdbool.h>

#include<openssl/err.h>
#include<openssl/evp.h>
#include<openssl/aes.h>
#include<openssl/rand.h>


typedef struct array
{
	int arr[256];
}crystal_poly;

typedef struct array1
{
	long long int arr[256];
}crystal_poly1;

//dilithium2
//int d=13,q=8380417,k=4,l=4,gamma1=131072,gamma2=95232,eta=2,beta=78,lamda=128,tau=39,omega=80;

//dilithium3
//int d=13,q=8380417,k=6,l=5,gamma1=524288,gamma2=261888,eta=4,beta=196,lamda=192,tau=49,omega=55;

//dilithium5
int d=13,q=8380417,k=8,l=7,gamma1=524288,gamma2=261888,eta=2,beta=120,lamda=256,tau=60,omega=75;


uint16_t messg_size=33;

int bitLength(int num);

void NTT(int f1[],long long f1_NTT[]);
void INTT(long long int f2[],int q,int tau);
void multiply_NTT(int a[],long long b[],long long out[]);

void shake_128(unsigned char rho[],unsigned char md_value[],unsigned int md_len,size_t size);
void shake_256(unsigned char input[],unsigned char md_value[],unsigned int md_len,size_t size);

// int coeff_from_three_bytes(int b0,int b1,int b2);
int coeff_from_three_bytes(unsigned char b0, unsigned char b1, unsigned char b2);

int coeff_from_half_bytes(int b);

void RejNTTPoly(unsigned char rho[],int j,int i,int out[]);
void RejBoundPoly(unsigned char rho_dash[],int out[]);

void power2Round(int t[],int t0[],int t1[],int d);

void int_to_bits(int input,int alpha,uint16_t out[]);
int bits_to_int(uint16_t input[],size_t size);

void bits_to_bytes(uint16_t input[],uint16_t out[],int c);    
void bytes_to_bits(uint16_t input[],uint16_t out[],size_t size);

void bit_pack(int w[],int a,int b,uint16_t out[]);
void simple_bit_pack(int input[],int b,uint16_t out[]);

void pk_Encode(uint16_t rho[],crystal_poly t1[],uint16_t pk[]);
void sk_encode(uint16_t rho[],unsigned char k1[],unsigned char tr[],crystal_poly s1[],crystal_poly s2[],crystal_poly t0[],uint16_t sk[]);

void keyGeneration(unsigned char seed[],uint16_t pk[],uint16_t sk[]);

//signature

void bit_unpack(uint16_t v[],int out[],int a,int b,size_t size);
void skDecode(uint16_t sk[],uint16_t rho[],uint16_t k1[],uint16_t tr[],crystal_poly s1[],crystal_poly s2[],crystal_poly t0[],uint16_t eta,uint16_t k,uint16_t l);

void ExpandMask(unsigned char rho_dash[],uint16_t neu,crystal_poly y[]);

void decompose(int r,int r2[]);
void High_bits(int input[],int out[]);

void w1_Encode(crystal_poly w1[],uint16_t bits[]);

void sample_in_ball(unsigned char c1_bar[],int c_out[],size_t size);

void Low_bits(int input[],int out[]);

int inf_norm_poly_ct0(long long int z[]);
int inf_norm_poly(int z[]);
int inf_norm_Grp(int z);

int make_hint(crystal_poly1 z[],crystal_poly r[],crystal_poly h[]);

void hint_bit_pack(crystal_poly h[],int y[]);
void sig_encode(unsigned char c_bar[],crystal_poly z[],crystal_poly h[],uint16_t sign[]);

void signature(uint16_t sk[],uint16_t message[],uint16_t sign[]);

//verification

void simple_bit_unpack(uint16_t v[],int out[],int b,size_t size);
void pk_decode(uint16_t pk[],uint16_t rho[],crystal_poly t1[]);

int hint_bit_unpack(uint16_t y[],crystal_poly h[]);
void sig_decode(uint16_t sign[],unsigned char c_bar[],crystal_poly z[],crystal_poly h[]);

void use_hint(int h[],long long r[],int w1[]);

bool verification(uint16_t pk[],uint16_t message[],uint16_t sign[]);

int main(){
    unsigned char seed[34];//={0xDF, 0xCC, 0x13, 0xCE, 0xD6, 0x97, 0x1E, 0xB1, 0xBF, 0x32, 0x43, 0xCB, 0x8E, 0xE8, 0x83, 0xFE, 0xA9, 0x67, 0x7D, 0x1E, 0x5D, 0xA8, 0xF3, 0x04, 0x6C, 0xFA, 0x43, 0x05, 0xDF, 0xB7, 0x91, 0x27,k,l};
    

    // if(k==4 && l==4)
    // {
    //     seed[32]="0x04";
    //     seed[33]="0x04";
    // }
    // else if(k==6 && l==5)
    // {
    //     seed[32]='0x06';
    //     seed[33]='0x05';
    // }
    // else if(k==8 && l==7)
    // {
    //     seed[32]='0x08';
    //     seed[33]='0x07';
    // }
    
    
    RAND_bytes(seed,32);
    seed[32]=k;
    seed[33]=l;

    uint16_t message[]={67, 11, 31, 70, 232, 125, 222, 154, 61, 5, 90, 125, 77, 106, 177, 39, 123, 45, 166, 237, 166, 66, 137, 100, 18, 18, 99, 145, 170, 43, 41, 175, 216};
    //uint16_t message_bits[4922*8];

    //size_t size=sizeof(message);
    //bytes_to_bits(message,message_bits,size);


    int len=bitLength(2*eta);
	int len1=bitLength(q-1);
    int len2=bitLength(gamma1-1);

    //32+32*k*(len1-d)
    uint16_t pk[32+32*k*(len1-d)];//={};
    
    //32+l*32*(1+len2)+omega+k
    uint16_t sign[(lamda/4)+l*32*(1+len2)+omega+k];//={};

    
    //128+32*((l+k)*len+d*k)
    uint16_t sk[128+32*((l+k)*len+d*k)];//={};
    bool ans=false;

    keyGeneration(seed,pk,sk);
    signature(sk,message,sign);
    ans=verification(pk,message,sign);

    if(ans==true)
    {
        printf("\nSignature verified\n");
    }
    else
    {
        printf("\nSignature is invalid\n");
    }

    return 0;
}

void keyGeneration(unsigned char seed[],uint16_t pk[],uint16_t sk[])
{
    printf("\nKey Generation\n");
    int len=bitLength(2*eta);

	int len1=bitLength(q-1);

    unsigned char shake256_out[128];

	size_t size=34;
	shake_256(seed,shake256_out,128,size); // Input: seed(34), Output: shake256_out(128)

    unsigned char rho[32],rho_dash[64],k1[32];

    for(int i=0;i<32;i++)
  	{
  		rho[i]=shake256_out[i];
  	} 	  	
  	// printf("\nRho :");
  	// for(int i=0;i<32;i++)
  	// {
  	// 	printf("%02x",rho[i]);
  	// }
	// printf("\n");

    int p=0;
	for(int i=32;i<96;i++)
	{
		rho_dash[p]=shake256_out[i];
		p++;
	}
	// printf("\nRho dash:");
  	// for(int i=0;i<64;i++)
  	// {
  	// 	printf("%02x",rho_dash[i]);
  	// }
	// printf("\n");

    p=0;
  	for(int i=96;i<128;i++)
  	{
  		k1[p]=shake256_out[i];
  		p++;
  	}
  	// printf("\nk :");
  	// for(int i=0;i<32;i++)
  	// {
  	// 	printf("%02x",k1[i]);
  	// }
  	// printf("\n");

    unsigned char rho1[34];
    uint16_t rho_uint[32];
	for(int i=0;i<32;i++)
  	{
  		rho1[i]=rho[i];
		rho_uint[i]=rho[i];
  	}

    //Expand A
	crystal_poly a[k][l];
		
	for(int i=0;i<k;i++)
	{
		for(int j=0;j<l;j++)
		{
			RejNTTPoly(rho1,j,i,a[i][j].arr);                            //KeyGeneration Algorithm line no. 3
		}
	}

    // printf("\nThe values of matrix A :\n");
	// for(int i=0;i<k;i++)
	// {
	// 	for(int j=0;j<l;j++)
	// 	{
	// 	    for(int var=0;var<256;var++)
	// 	    {
	// 		    printf("%d ",a[i][j].arr[var]);
	// 	    }
	// 		printf("\n");
	// 	}
	// 	printf("\n");
	// }


	//Expand S
	crystal_poly s1[l],s2[k];
	crystal_poly1 s1_NTT[l];

    unsigned char rho_dash1[66];
    for(int i=0;i<64;i++)
  	{
  		rho_dash1[i]=rho_dash[i];
  	}
		


    //Expand S    
	for(int r=0;r<l;r++)
	{
		rho_dash1[64]=r;
		rho_dash1[65]=0;
		RejBoundPoly(rho_dash1,s1[r].arr);                             //KeyGeneration Algorithm line no. 4
	}

	printf("\nThe values of s1 (Modified changed) :\n");
	for(int r=0;r<l;r++)
	{
		for(int var=0;var<256;var++)
		    {
			    printf("%d ",s1[r].arr[var]);
		    }
			printf("\n\n");
	}
	printf("\n");

    for(uint16_t r=0;r<k;r++)
	{
		rho_dash1[64]=r+l;
		rho_dash1[65]=0;
		RejBoundPoly(rho_dash1,s2[r].arr);
	}
	// printf("\nThe values of s2 :\n");
	// for(uint16_t r=0;r<k;r++)
	// {
	// 	for(int var=0;var<256;var++)
	// 	    {
	// 		    printf("%d ",s2[r].arr[var]);
	// 	    }
	// 		printf("\n\n");
	// }
	// printf("\n");
	
	for (int i = 0; i < l; i++)
	{
		NTT(s1[i].arr,s1_NTT[i].arr);
	}
	// printf("\ns1 NTT:\n");
	// for(uint16_t r=0;r<l;r++)
	// {
	// 	for(int var=0;var<256;var++)
	// 	    {
	// 		    printf("%lld ",s1_NTT[r].arr[var]);
	// 	    }
	// 		printf("\n");
	// }
	// printf("\n");


    crystal_poly t[k];
	crystal_poly1 result[k];

    //Multiplication of A and S1
	long long mul[256];
	for(int i=0;i<k;i++)
	{
		for(int p=0;p<256;p++)
		{
			result[i].arr[p]=0;
		}
		for(int j=0;j<l;j++)
		{
			multiply_NTT(a[i][j].arr,s1_NTT[j].arr,mul);
			for(int m=0;m<256;m++)
			{
				result[i].arr[m]+=mul[m];
				//result[i].arr[m]=result[i].arr[m] % q;
				//printf("%d ",result[i].arr[m]);

			}
		}
	}
	
	for (int i = 0; i < k; i++)
	{
		INTT(result[i].arr,q,1753);
	}

    for(int i=0;i<k;i++)
	{
		for(int j=0;j<256;j++)
		{
			t[i].arr[j]=result[i].arr[j] + s2[i].arr[j];                   //KeyGeneration Algorithm line no. 5     
		}
	}

	// printf("\nThe values of t :\n");
	// for(uint16_t r=0;r<k;r++)
	// {
	// 	for(int var=0;var<256;var++)
	// 	    {
	// 		    printf("%d ",t[r].arr[var]);
	// 	    }
	// 		printf("\n\n");
	// }
	// printf("\n");
	
	crystal_poly t0[k],t1[k];

	for(int i=0;i<k;i++)
	{
		power2Round(t[i].arr,t0[i].arr,t1[i].arr,d);                     //KeyGeneration Algorithm line no. 6
	}
	
	// printf("\nThe values of t0 :\n");
	// for(uint16_t r=0;r<k;r++)
	// {
	// 	for(int var=0;var<256;var++)
	// 	    {
	// 		    printf("%d ",t0[r].arr[var]);
	// 	    }
	// 		printf("\n\n");
	// }
	// printf("\n");

    // printf("\nThe values of t1 :\n");
	// for(uint16_t r=0;r<k;r++)
	// {
	// 	for(int var=0;var<256;var++)
	// 	    {
	// 		    printf("%d ",t1[r].arr[var]);
	// 	    }
	// 		printf("\n\n");
	// }
	// printf("\n");

	//pk Encode
	pk_Encode(rho_uint,t1,pk);                              //KeyGeneration Algorithm line no. 8

	// printf("\nPk :\n");
	// for (int i = 0; i < 32+32*k*(23-d); i++)
	// {
	// 	printf("%02x",pk[i]);
	// }

	unsigned char pk1[32+32*k*(23-d)];
	for (int i = 0; i < 32+32*k*(23-d); i++)
	{
		pk1[i]=pk[i];
	}

    unsigned char tr[64];

    size=sizeof(pk1);	
	shake_256(pk1,tr,64,size);                      //KeyGeneration Algorithm line no. 9

	// printf("\ntr :");
    // for(int i=0;i<64;i++)
    // {
    //     printf("%02x ",tr[i]);
    // }
    // printf("\n");

	// printf("\nk1 :");
    // for(int i=0;i<32;i++)
    // {
    //     printf("%02x",k1[i]);
    // }
    // printf("\n");

	// printf("\nThe values of s1 :\n");
	// for(int r=0;r<l;r++)
	// {
	// 	for(int var=0;var<256;var++)
	// 	    {
	// 		    printf("%d ",s1[r].arr[var]);
	// 	    }
	// 		printf("\n\n");
	// }
	// printf("\n");

	// printf("\nThe values of s2 :\n");
	// for(uint16_t r=0;r<k;r++)
	// {
	// 	for(int var=0;var<256;var++)
	// 	    {
	// 		    printf("%d ",s2[r].arr[var]);
	// 	    }
	// 		printf("\n\n");
	// }
	// printf("\n");

	sk_encode(rho_uint,k1,tr,s1,s2,t0,sk);                     //KeyGeneration Algorithm line no. 10
		
	// printf("\nSk :\n");
	// for (int i = 0; i < 128+32*((l+k)*len+d*k); i++)
	// {
	// 	printf("%02x",sk[i]);
	// }    
	// printf("\n");

}

void signature(uint16_t sk[],uint16_t message[],uint16_t sign[])
{
    printf("\n Signing \n");
    int len1=bitLength(gamma1-1);

    uint16_t rho[32],k1[32],tr[64];
    crystal_poly s1[l],s2[k],t0[k];
    crystal_poly1 s1_NTT[l],s2_NTT[k],t0_NTT[k];

    skDecode(sk,rho,k1,tr,s1,s2,t0,eta,k,l);                   //Signature Algorithm line no. 1

    // printf("\nRho :");
    // for(int i=0;i<32;i++)
    // {
    //     printf("%02x",rho[i]);
    // }

    // printf("\nk1 :");
    // for(int i=0;i<32;i++)
    // {
    //     printf("%02x",k1[i]);
    // }
    // printf("\n");

    // printf("\ntr :");
    // for(int i=0;i<64;i++)
    // {
    //     printf("%02x",tr[i]);
    // }
    // printf("\n");

    // printf("\ns1 :\n");
    // for(int i=0;i<l;i++)
    // {
    //     for(int j=0;j<256;j++)
    //     {
    //         printf("%d ",s1[i].arr[j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    // printf("\ns2 :\n");
    // for(int i=0;i<k;i++)
    // {
    //     for(int j=0;j<256;j++)
    //     {
    //         printf("%d ",s2[i].arr[j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    // printf("\nt0 :\n");
    // for(int i=0;i<k;i++)
    // {
    //     for(int j=0;j<256;j++)
    //     {
    //         printf("%d ",t0[i].arr[j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");


    for(int i=0;i<l;i++)
    {
        NTT(s1[i].arr,s1_NTT[i].arr);                   //Signature Algorithm line no. 2
    }

    for(int i=0;i<k;i++)
    {
        NTT(s2[i].arr,s2_NTT[i].arr);                   //Signature Algorithm line no. 3
    }
    for(int i=0;i<k;i++)
    {
        NTT(t0[i].arr,t0_NTT[i].arr);                   //Signature Algorithm line no. 4
    }

    //Expand A
    unsigned char rho1[34];
    for (int i = 0; i < 32; i++)
    {
        rho1[i]=rho[i];
    }
    
	crystal_poly a[k][l];
		
	for(int i=0;i<k;i++)
	{
		for(int j=0;j<l;j++)
		{
			RejNTTPoly(rho1,j,i,a[i][j].arr);           //Signature Algorithm line no. 5
		}
	}

    // printf("\nThe values of matrix A :\n");
	// for(int i=0;i<k;i++)
	// {
	// 	for(int j=0;j<l;j++)
	// 	{
	// 	    for(int var=0;var<256;var++)
	// 	    {
	// 		    printf("%d ",a[i][j].arr[var]);
	// 	    }
	// 		printf("\n");
	// 	}
	// 	printf("\n");
	// }


    unsigned char shake_256ip[64+(messg_size)],neu[64];
    int count=0;
    for (int i = 0; i < 64; i++)
    {
        shake_256ip[i]=tr[i];
        count++;
    }
    for (int i = 0; i < messg_size; i++)
    {
        shake_256ip[count]=message[i];
        count++;
    }

    size_t size=sizeof(shake_256ip);
    shake_256(shake_256ip,neu,64,size);             //Signature Algorithm line no. 6

    // printf("\nneu :\n");
    // for (int i = 0; i < 64; i++)
    // {
    //     printf("%02x",neu[i]);
    // }
    // printf("\n");

    // 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    unsigned char rnd[32]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};              //For Determinsitic Varient
	//RAND_bytes(rnd,32);

	unsigned char shake_256ip1[128],rho_dash[64];
    count=0;
    for(int i=0;i<32;i++)
    {
        shake_256ip1[i]=k1[i];
        count++;
    }
    for(int i=0;i<32;i++)
    {
        shake_256ip1[count]=rnd[i];
        count++;
    }
    for(int i=0;i<64;i++)
    {
        shake_256ip1[count]=neu[i];
        count++;
    }

    size=sizeof(shake_256ip1);
    shake_256(shake_256ip1,rho_dash,64,size);       //Signature Algorithm line no. 7

    // printf("\nrho dash :\n");
    // for (int i = 0; i < 64; i++)
    // {
    //     printf("%02x",rho_dash[i]);
    // }
    // printf("\n");

    uint16_t counter_k=0;

    //int z=-1,h=-1;
    crystal_poly z[l],h[k];
    for (int i = 0; i < l; i++)
    {
        for (int j = 0; j < 256; j++)
        {
            z[i].arr[j]=-1;                        //Initialize z to blank symbol
        }
    }
    
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < 256; j++)
        {
            h[i].arr[j]=-1;                       //Initialize h to blank symbol
        }
    }

    int count_z=0,count_h=0;

    unsigned char c_bar[(2*lamda)/8];

    crystal_poly w[k],w1[k];

    crystal_poly y[l];
    crystal_poly1 y_NTT[l];
    crystal_poly1 cs1[l],cs2[k],ct0[k];


    int no_of_ones=0;
    while (count_z == 0 && count_h == 0){

        ExpandMask(rho_dash,counter_k,y);       //Signature Algorithm line no. 11
       

        // printf("\ny =\n");
        // for (int i = 0; i < l; i++)
        // {
        //     printf("\n");
        //     for (int j = 0; j < 256; j++)
        //     {
        //         printf("%d ",y[i].arr[j]);
        //     }
        //     printf("\n");
        // }
    	// 	printf("\n");

        

        for(int i=0;i<l;i++)
        {
            NTT(y[i].arr,y_NTT[i].arr);
        }

        
        //w = A * y
        crystal_poly1 result1[k];
        
        //Multiplication of A and y
	    long long mul[256];
	    for(int i=0;i<k;i++)
	    {
		    for(int p=0;p<256;p++)
		    {
			    result1[i].arr[p]=0;
		    }
		    for(int j=0;j<l;j++)
		    {
			    multiply_NTT(a[i][j].arr,y_NTT[j].arr,mul);
			    for(int m=0;m<256;m++)
			    {
				    result1[i].arr[m]+=mul[m];
				//result[i].arr[m]=result[i].arr[m] % q;
				//printf("%d ",result[i].arr[m]);

			    }
		    }
	    }

        for(int i=0;i<k;i++)
        {
           INTT(result1[i].arr,q,1753); 
        }


        for (int i = 0; i < k; i++)
        {
            for (int j = 0; j < 256; j++)
            {
                w[i].arr[j]=result1[i].arr[j];                  //Signature Algorithm line no. 12
            } 
        }

        // printf("\nw =\n");
        // for (int i = 0; i < k; i++)
        // {
        //     for (int j = 0; j < 256; j++)
        //     {
        //         printf("%d ",w[i].arr[j]);
        //     }
        //     printf("\n");
        // }
    	// 	printf("\n");


        for(int i=0;i<k;i++)
        {
           High_bits(w[i].arr,w1[i].arr);                       //Signature Algorithm line no. 13
        }
        
        // printf("\nw1 =\n");
        // for (int i = 0; i < k; i++)
        // {
        //     for (int j = 0; j < 256; j++)
        //     {
        //         printf("%d ",w1[i].arr[j]);
        //     }
        //     printf("\n");
        // }
    	// 	printf("\n");


        

        int x=(q-1)/(2*gamma2)- 1 ;
        int len2=bitLength(x);        
        uint16_t w1_Encode_out[32*k*len2];
        w1_Encode(w1,w1_Encode_out); 

        // printf("\nw1_encode_out =\n");
        // for (int i = 0; i < 32*k*len2; i++)
        // {
            
        //         printf("%02x",w1_Encode_out[i]);

        // }
    	// printf("\n");

        
        unsigned char shake_256ip_2[64+32*k*len2];
        size=sizeof(shake_256ip_2);
        count=0;
        for (int i = 0; i < 64; i++)
        {
            shake_256ip_2[count]=neu[i];
            count++;
        }
        for (int i = 0; i < 32*k*len2; i++)
        {
            shake_256ip_2[count]=w1_Encode_out[i];
            count++;
        }

        shake_256(shake_256ip_2,c_bar,(2*lamda)/8,size);                    //Signature Algorithm line no. 15
        
        unsigned char c1_bar[32];//c2_bar[];
        for (int i = 0; i < 32; i++)
        {
           c1_bar[i]=c_bar[i];
        }

        int c[256];

        size=sizeof(c_bar);
        sample_in_ball(c_bar,c,size);                                   //Signature Algorithm line no. 16

        long long c_NTT[256];
        NTT(c,c_NTT);                                                  //Signature Algorithm line no. 17

        // printf("\nc NTT =\n");
        // for (int i = 0; i < 256; i++)
        // {  
        //     printf("%d ",c_NTT[i]);
        // }
    	// printf("\n");

        // printf("\ns1 NTT :\n");
        // for(int i=0;i<l;i++)
        // {
        //     for(int j=0;j<256;j++)
        //     {
        //         printf("%d ",s1_NTT[i].arr[j]);
        //     }
        //     printf("\n");
        // }
        // printf("\n");
        
        for (int i = 0; i < l; i++)
        {
            for (int j = 0; j < 256; j++)
            {
                cs1[i].arr[j]=(c_NTT[j]*s1_NTT[i].arr[j])%q;
            }
        }
        
        for (int i = 0; i < l; i++)
        {
            INTT(cs1[i].arr,q,1753);                                    //Signature Algorithm line no. 18
        }

        // printf("\ncs1 =\n");
        // for (int i = 0; i < l; i++)
        // {
        //     for (int j = 0; j < 256; j++)
        //     {
        //         printf("%d ",cs1[i].arr[j]);
        //     }
        //     printf("\n");
        // }
    	// 	printf("\n");
        
        for (int i = 0; i < k; i++)
        {
            for (int j = 0; j < 256; j++)
            {
                cs2[i].arr[j]=(c_NTT[j]*s2_NTT[i].arr[j])%q;
            }
        }
        
        for (int i = 0; i < k; i++)
        {
            INTT(cs2[i].arr,q,1753);                                  //Signature Algorithm line no. 19
        }

        // printf("\ncs2 =\n");
        // for (int i = 0; i < k; i++)
        // {
        //     for (int j = 0; j < 256; j++)
        //     {
        //         printf("%d ",cs2[i].arr[j]);
        //     }
        //     printf("\n");
        // }
    	// printf("\n");

        for (int i = 0; i < l; i++)
        {
            for (int j = 0; j < 256; j++)
            {
                z[i].arr[j] = (y[i].arr[j] + cs1[i].arr[j])%q;                  //Signature Algorithm line no. 20
                if(z[i].arr[j] < 0)
                {
                    z[i].arr[j]=z[i].arr[j]+q;
                }
            }           
        }

        // printf("\nz =\n");
        // for (int i = 0; i < l; i++)
        // {
        //     for (int j = 0; j < 256; j++)
        //     {
        //         printf("%d ",z[i].arr[j]);
        //     }
        //     printf("\n");
        // }
    	// printf("\n");
        
        crystal_poly low_bits[k],r0[k];
        for (int i = 0; i < k; i++)
        {
            for (int j = 0; j < 256; j++)
            {
                low_bits[i].arr[j] = w[i].arr[j] - cs2[i].arr[j];
            } 
        }

        for (int i = 0; i < k; i++)
        {
            Low_bits(low_bits[i].arr,r0[i].arr);                             //Signature Algorithm line no. 21
        }

        // printf("\nr0 =\n");
        // for (int i = 0; i < k; i++)
        // {
        //     for (int j = 0; j < 256; j++)
        //     {
        //         printf("%d ",r0[i].arr[j]);
        //     }
        //     printf("\n");
        // }
    	// 	printf("\n");

        
        // for (int i = 0; i < k; i++)
        // {
        //     for (int j = 0; j < 256; j++)
        //     {
        //         r0[i].arr[j]=r0[i].arr[j] % q;
        //         if(r0[i].arr[j]<0)
        //         {
        //             r0[i].arr[j]+=q;
        //         }
        //     }
        // }


        int var1=gamma1-beta;
        int var2=gamma2-beta;

        int z_mod[l];
        for (int i = 0; i < l; i++)
        {
            z_mod[i]=inf_norm_poly(z[i].arr);                           //centered modular reduction operation
        }
        

        int z_mod_max=z_mod[0];
        for (int i = 1; i < l; i++)
        {
            if(z_mod[i] > z_mod_max)
            {
                z_mod_max=z_mod[i];                                     //maximum of z mod q
            }
        }

        int r0_mod[k];
        for (int i = 0; i < k; i++)
        {
            r0_mod[i]=inf_norm_poly(r0[i].arr);                         //centered modular reduction operation
        }

        int r0_mod_max=r0_mod[0];
        for (int i = 1; i < k; i++)
        {
            if(r0_mod[i] > r0_mod_max)
            {
                r0_mod_max=r0_mod[i];                                   //maximum of r mod q
            }
        }

        //printf("\nz_mod_max=%d\n",z_mod_max);
        //printf("\nr0_mod_max=%d\n",r0_mod_max);

        //break;
       
        if(z_mod_max >= var1 || r0_mod_max >= var2)
        {
            for (int i = 0; i < l; i++)
            {
                for (int j = 0; j < 256; j++)
                {
                    z[i].arr[j]=-1;
                }
            }
    
            for (int i = 0; i < k; i++)
            {
                for (int j = 0; j < 256; j++)
                {
                    h[i].arr[j]=-1;
                }
            }
        }        
        else
        {
            for (int i = 0; i < k; i++)
            {
                for (int j = 0; j < 256; j++)
                {
                    ct0[i].arr[j]=c_NTT[j]*t0_NTT[i].arr[j];                    
                }
            }

            for (int i = 0; i < k; i++)
            {
                INTT(ct0[i].arr,q,1753);                               //Signature Algorithm line no. 25
            }

            // printf("\nct0 =\n");
            // for (int i = 0; i < k; i++)
            // {
            //     for (int j = 0; j < 256; j++)
            //     {
            //         printf("%d ",ct0[i].arr[j]);
            //     }
            //     printf("\n");
            // }
    		// printf("\n");


            crystal_poly1 ct0_neg[k];
            for (int i = 0; i < k; i++)
            {
                for (int j = 0; j < 256; j++)
                {
                    ct0_neg[i].arr[j]=q-ct0[i].arr[j];
                }
            }

            // printf("\nct negative =\n");
            // for (int i = 0; i < k; i++)
            // {
            //     for (int j = 0; j < 256; j++)
            //     {
            //         printf("%d ",ct0_neg[i].arr[j]);
            //     }
            //     printf("\n");
            // }
    		// printf("\n");

        //     printf("\nw =\n");
        // for (int i = 0; i < k; i++)
        // {
        //     for (int j = 0; j < 256; j++)
        //     {
        //         printf("%d ",w[i].arr[j]);
        //     }
        //     printf("\n");
        // }
    	// 	printf("\n");



        //     printf("\ncs2 =\n");
        //     for (int i = 0; i < k; i++)
        //     {
        //         for (int j = 0; j < 256; j++)
        //         {
        //             printf("%d ",cs2[i].arr[j]);
        //         }
        //         printf("\n");
        //     }
    	//     printf("\n");

        //     printf("\nct0 =\n");
        //     for (int i = 0; i < k; i++)
        //     {
        //         for (int j = 0; j < 256; j++)
        //         {
        //             printf("%d ",ct0[i].arr[j]);
        //         }
        //         printf("\n");
        //     }
    	// 	printf("\n");

            crystal_poly result[k];
            for (int i = 0; i < k; i++)
            {
                for (int j = 0; j < 256; j++)
                {
                    result[i].arr[j] =( w[i].arr[j] - cs2[i].arr[j] + ct0[i].arr[j] ) % q;
                }
            }
            
            // printf("\nw-cs2+ct0 =\n");
            // for (int i = 0; i < k; i++)
            // {
            //     for (int j = 0; j < 256; j++)
            //     {
            //         printf("%d ",result[i].arr[j]);
            //     }
            //     printf("\n");
            // }
    		// printf("\n");

            // printf("\nh beforemake hint=\n");
            // for (int i = 0; i < k; i++)
            // {
            //     for (int j = 0; j < 256; j++)
            //     {
            //         printf("%d ",h[i].arr[j]);
            //     }
            //     printf("\n");
            // }
    		// printf("\n");


            //printf("\nno of ones before make hint=%d\n",no_of_ones);
            no_of_ones=make_hint(ct0_neg,result,h);                                     //Signature Algorithm line no. 26
            //printf("\nno of ones after make hint=%d\n",no_of_ones);

            int ct0_mod[k];
            for (int i = 0; i < k; i++)
            {
                ct0_mod[i]=inf_norm_poly_ct0(ct0[i].arr);                               //centered modular reduction operation
            }
            int ct0_mod_max=ct0_mod[0];
            for (int i = 1; i < k; i++)
            {
                if(ct0_mod[i] > ct0_mod_max)
                {
                    ct0_mod_max=ct0_mod[i];                                             //maximum of ct0 mod q
                }
            }

            //printf("\nct0_mod_max=%d\n",ct0_mod_max);

            if(ct0_mod_max >= gamma2 || no_of_ones > omega)
            {
                for (int i = 0; i < l; i++)
                {
                    for (int j = 0; j < 256; j++)
                    {
                        z[i].arr[j]=-1;
                    }
                }
    
                for (int i = 0; i < k; i++)
                {
                    for (int j = 0; j < 256; j++)
                    {
                        h[i].arr[j]=-1;
                    }
                }
            }

           // break;

        }

        counter_k = counter_k + l;                      //Signature Algorithm line no. 31
        for (int i = 0; i < l; i++)
        {
            for (int j = 0; j < 256; j++)
            {
                if(z[i].arr[j] != -1)
                {
                    count_z=1;                          //loop exit condition
                    break;
                }
            }
        }

        for (int i = 0; i < k; i++)
        {
            for (int j = 0; j < 256; j++)
            {
                if(h[i].arr[j] != -1)
                {
                    count_h=1;                          //loop exit condition
                    break;
                }
            }
        }   

        if (z_mod_max >= var1) {
            printf("❌ Rejected: z_mod_max (%d) >= gamma1 - beta (%d)\n", z_mod_max, var1);
        }
        else{
            printf("✅ Accepted: z_mod_max (%d) < gamma1 - beta (%d)\n", z_mod_max, var1);
        }
        if (r0_mod_max >= var2) {
            printf("❌ Rejected: r0_mod_max (%d) >= gamma2 - beta (%d)\n", r0_mod_max, var2);
        }
        else{
            printf("✅ Accepted: r0_mod_max (%d) < gamma2 - beta (%d)\n", r0_mod_max, var2);
        }
        
    
    }
 
    for (int i = 0; i < l; i++)
    {
        for (int j = 0; j < 256; j++)
        {
            z[i].arr[j]=inf_norm_Grp(z[i].arr[j]);                      //centered modular reduction operation
        }
    }
    
    sig_encode(c_bar,z,h,sign);                         //Signature Algorithm line no. 33

    // printf("\nsign :\n");
    // for (int i = 0;i<(lamda/4)+l*32*(1+len1)+omega+k; i++)
    // {
    //     printf("%02x",sign[i]);
    // }
    printf("\n");

    // Debugging: Check where rejection happens



}

bool verification(uint16_t pk[],uint16_t message[],uint16_t sign[]){
    printf("\nVerification\n");
    int len=bitLength(q-1);
    int len1=bitLength(gamma1-1);

    uint16_t rho[32];
    crystal_poly t1[k];

    pk_decode(pk,rho,t1);                               //Verification Algorithm line no. 1

    // printf("\nRho :");
    // for(int i=0;i<32;i++)
    // {
    //     printf("%02x",rho[i]);
    // }
	// printf("\n");


    // printf("\nThe values of t1 :\n");
	// for(int r=0;r<k;r++)
	// {
	// 	for(int var=0;var<256;var++)
	// 	    {
	// 		    printf("%d ",t1[r].arr[var]);
	// 	    }
	// 		printf("\n\n");
	// }
	// printf("\n");

    unsigned char c_bar[(2*lamda)/8];                //c tilde

    crystal_poly z[l],h[k];

    sig_decode(sign,c_bar,z,h);                      //Verification Algorithm line no. 2

    printf("\nc_bar :");
    for(int i=0;i<(2*lamda)/8;i++)
    {
        printf("%02x",c_bar[i]);
    }
	printf("\n");

    // printf("\nz :\n");
	// for(int r=0;r<l;r++)
	// {
	// 	for(int var=0;var<256;var++)
	// 	    {
	// 		    printf("%d ",z[r].arr[var]);
	// 	    }
	// 		printf("\n\n");
	// }
	// printf("\n");

    // printf("\nh :\n");
	// for(int r=0;r<k;r++)
	// {
	// 	for(int var=0;var<256;var++)
	// 	    {
	// 		    printf("%d, ",h[r].arr[var]);
	// 	    }
	// 		printf("\n\n");
	// }
	// printf("\n");


    int count=0;
    for (int i = 0; i < k; i++)
{
        for (int j = 0; j < 256; j++)
        {
            if(h[i].arr[j]== -1)
            {
                count=1;
                break;                                      
            }
        }
    }

    if(count == 1)
    {
        return false;
    }
    else
    {
        crystal_poly a[k][l];
		
        unsigned char rho1[34];
        for (int i = 0; i < 32; i++)
        {
            rho1[i]=rho[i];
        }

	    for(int i=0;i<k;i++)
	    {
		    for(int j=0;j<l;j++)
		    {
			    RejNTTPoly(rho1,j,i,a[i][j].arr);                           //Verification Algorithm line no. 5
		    }
	    }

        // printf("\nThe values of matrix A :\n");
	    // for(int i=0;i<k;i++)
	    // {
		//     for(int j=0;j<l;j++)
		//     {
		//         for(int var=0;var<256;var++)
		//         {
		// 	        printf("%d ",a[i][j].arr[var]);
		//         }
		// 	    printf("\n");
		//     }
		//     printf("\n");
	    // }

        unsigned char tr[64];
        unsigned char pk1[(32+32*k*(23-d))];

        // uint16_t pk_copy[32+32*k*(23-d)],pk_bits[(32+32*k*(23-d))*8];

        // for (int i = 0; i < 32+32*k*(23-d); i++)
	    // {
		//     pk_copy[i]=pk[i];
	    // }


        // size_t size1=sizeof(pk_copy);


        // bytes_to_bits(pk_copy,pk_bits,size1);


	    for(int i = 0; i < (32+32*k*(23-d)); i++)
	    {
		    pk1[i]=pk[i];
	    }
        size_t size=sizeof(pk1);
	    shake_256(pk1,tr,64,size);                                      //Verification Algorithm line no. 6

        // printf("\ntr :");
        // for(int i=0;i<64;i++)
        // {
        //     printf("%d ",tr[i]);
        // }
	    // printf("\n");

        // uint16_t tr_int[64],tr_bits[512];
        // size=sizeof(tr_int);
        // for (int i = 0; i < 64; i++)
        // {
        //     tr_int[i]=tr[i];
        // }
        
        // bytes_to_bits(tr_int,tr_bits,size);

        // printf("\ntr bits = \n");
        // for (int i = 0; i < 512; i++)
        // {
        //     printf("%d ",tr_bits[i]);
            
        // }
        // printf("\n");
        
        //printf("\nmessage size = %d \n",messg_size);

        unsigned char shake_256ip[64+messg_size],neu[64];
        int count=0;
        for (int i = 0; i < 64; i++)
        {
            shake_256ip[i]=tr[i];
            count++;
        }
        for (int i = 0; i < messg_size; i++)
        {
            shake_256ip[count]=message[i];
            count++;
        }

        // printf("\nshake256 input = \n");
        // for (int i = 0; i < 64 + messg_size; i++)
        // {
        //     printf("%d ",shake_256ip[i]);
            
        // }
        // printf("\n");


        size=sizeof(shake_256ip);
        shake_256(shake_256ip,neu,64,size);                               //Verification Algorithm line no. 7

        // printf("\nneu :");
        // for(int i=0;i<64;i++)
        // {
        //     printf("%02x",neu[i]);
        // }
	    // printf("\n");

        int c[256];
        unsigned char c1_bar[32];
        for(int i=0;i<32;i++)
        {
            c1_bar[i]=c_bar[i];
        }

        // printf("\nc1 bar :");
        // for(int i=0;i<32;i++)
        // {
        //     printf("0x%02x, ",c1_bar[i]);
        // }
	    // printf("\n");

        size=sizeof(c_bar);
        sample_in_ball(c_bar,c,size);                                     //Verification Algorithm line no. 8

        // printf("\nc :");
        // for(int i=0;i<256;i++)
        // {
        //     printf("[%d] -> %d ,",i,c[i]);
        // }
	    // printf("\n");

        for (int i = 0; i < k; i++)
        {
            for (int j = 0; j < 256; j++)
            {
                t1[i].arr[j]=t1[i].arr[j] * pow(2,d);
            }
        }

        crystal_poly1 t1_NTT[k];
        for (int i = 0; i < k; i++)
        {
            NTT(t1[i].arr,t1_NTT[i].arr);
        }

        // printf("\nt1 NTT=\n");
        //     for (int i = 0; i < k; i++)
        //     {
        //         for (int j = 0; j < 256; j++)
        //         {
        //             printf("%d ",t1_NTT[i].arr[j]);
        //         }
        //         printf("\n");
        //     }
    	// 	printf("\n");   

        long long c_NTT[256];
        NTT(c,c_NTT);

        // printf("\nC NTT =\n");
        // for (int i = 0; i < 256; i++)
        // {
        //     printf("%d ",c_NTT[i]);
        // }
        // printf("\n");


        crystal_poly1 result1[k];
        for (int i = 0; i < k; i++)
        {
            for (int j = 0; j < 256; j++)
            {
                result1[i].arr[j]=(t1_NTT[i].arr[j] * c_NTT[j])%q;
            }
        }

        crystal_poly1 z_NTT[l];
        for (int i = 0; i < l; i++)
        {
            NTT(z[i].arr,z_NTT[i].arr);
        }

        // printf("\nz NTT=\n");
        //     for (int i = 0; i < l; i++)
        //     {
        //         for (int j = 0; j < 256; j++)
        //         {
        //             printf("%d ",z_NTT[i].arr[j]);
        //         }
        //         printf("\n");
        //     }
    	// 	printf("\n");

        crystal_poly1 result2[k];
        long long mul[256];
	    for(int i=0;i<k;i++)
	    {
		    for(int p=0;p<256;p++)
		    {
			    result2[i].arr[p]=0;
		    }
		    for(int j=0;j<l;j++)
		    {
			    multiply_NTT(a[i][j].arr,z_NTT[j].arr,mul);                     //multiplication of A matrix and z
			    for(int m=0;m<256;m++)
			    {
				    result2[i].arr[m]+=mul[m];
				result2[i].arr[m]=result2[i].arr[m] % q;
				//printf("%d ",result[i].arr[m]);

			    }
		    }
	    }

        // printf("\na*zntt=\n");
        //     for (int i = 0; i < k; i++)
        //     {
        //         for (int j = 0; j < 256; j++)
        //         {
        //             printf("%d ",result2[i].arr[j]);
        //         }
        //         printf("\n");
        //     }
    	// 	printf("\n");

        // printf("\nresult1=\n");
        //     for (int i = 0; i < k; i++)
        //     {
        //         for (int j = 0; j < 256; j++)
        //         {
        //             printf("%d ",result1[i].arr[j]);
        //         }
        //         printf("\n");
        //     }
    	// 	printf("\n");


         

        crystal_poly1 w_approx[k];
        for (int i = 0; i < k; i++)
        {
            for (int j = 0; j < 256; j++)
            {
                w_approx[i].arr[j]=result2[i].arr[j] - result1[i].arr[j];
            }
        }

        for (int i = 0; i < k; i++)
        {
            INTT(w_approx[i].arr,q,1753);                               //Verification Algorithm line no. 9
        }

        // printf("\nw approx=\n");
        //     for (int i = 0; i < k; i++)
        //     {
        //         for (int j = 0; j < 256; j++)
        //         {
        //             printf("%d ",w_approx[i].arr[j]);
        //         }
        //         printf("\n");
        //     }
    	// 	printf("\n");

        crystal_poly w1_dash[k];

        for (int i = 0; i < k; i++)
        {
            use_hint(h[i].arr,w_approx[i].arr,w1_dash[i].arr);              //Verification Algorithm line no. 10
        }

        // printf("\nw dash=\n");
        //     for (int i = 0; i < k; i++)
        //     {
        //         for (int j = 0; j < 256; j++)
        //         {
        //             printf("%d ",w1_dash[i].arr[j]);
        //         }
        //         printf("\n");
        //     }
    	// 	printf("\n");
        
        int x=(q-1)/(2*gamma2) - 1;
        int len2=bitLength(x);        
        uint16_t w1_Encode_out[32*k*len2];
        w1_Encode(w1_dash,w1_Encode_out);

        // printf("\nw1 encode out :\n");
        // for(int i=0;i<32*k*len2;i++)
        // {
        //     printf("%02x",w1_Encode_out[i]);
        // }
	    // printf("\n");


        //printf("\nafter w1 encode\n");
        unsigned char shake256_ip1[64+32*k*len2];
        count=0;
        for (int i = 0; i < 64; i++)
        {
            shake256_ip1[count]=neu[i];
            count++;
        }
        for (int i = 0; i < 32*k*len2; i++)
        {
            shake256_ip1[count]=w1_Encode_out[i];
            count++;
        }
        
        // uint16_t shake_ip_bytes[64+32*k*len2],shake_ip_bits[(64+32*k*len2)*8];
        // size=sizeof(shake_ip_bytes);
        // for (int i = 0; i < 64+32*k*len2; i++)
        // {
        //     shake_ip_bytes[i]=shake256_ip1[i];
        // }

        // bytes_to_bits(shake_ip_bytes,shake_ip_bits,size);



        unsigned char c_bar_dash[(2*lamda)/8];
        // for (int i = 0; i < (64+32*k*len2)*8; i++)
        // {
        //     shake_ip_bits1[i]=shake_ip_bits[i];
        // }


        size=sizeof(shake256_ip1);
        shake_256(shake256_ip1,c_bar_dash,(2*lamda)/8,size);                        //Verification Algorithm line no. 12

        printf("\nc_bar dash:");
        for(int i=0;i<(2*lamda)/8;i++)
        {
            printf("%02x",c_bar_dash[i]);
        }
	    printf("\n");


        int var1=gamma1-beta;
        int z_mod[l];
        for (int i = 0; i < l; i++)
        {
            z_mod[i]=inf_norm_poly(z[i].arr);                                     //centered modular reduction operation     
        }
        int z_mod_max=z_mod[0];
        for (int i = 1; i < k; i++)
        {
            if(z_mod[i] > z_mod_max)
            {
                z_mod_max=z_mod[i];                                               //maximum of z mod q
            }
        }

        count=0;
        for (int i = 0; i < (2*lamda)/8; i++)
        {
            if(c_bar_dash[i]==c_bar[i])
            {
                count++;
            }
        }

        // int count1=0;
        // for (int i = 0; i < k; i++)
        // {
        //     for (int j = 0; j < 256; j++)
        //     {
        //         if(h[i].arr[j] == 1)
        //         {
        //             count1++;
        //         }
        //     }
        // }

        if (z_mod_max>var1)
        {
            printf("\n z is too large \n");
        }        

        if( count!=(2*lamda)/8)
        {
            printf("\n modify signature \n");
        }
        

        if(z_mod_max<var1 && count==(2*lamda)/8 )                               //Verification Algorithm line no. 13
        {
            printf("\nno modification\n");
            return true;
        }
        else
        {
            return false;
        }

    }

}




//Primitives

int bitLength(int num)
{
    int length = 0;
    
    while (num != 0) {
        num = num >> 1;
        length++;
    }
    
    return length;
}


int bitreversal(int n) {
    int reversed = 0;
    for (int i = 0; i < 8; ++i) {
        if (n & (1 << i)) {
            reversed |= (1 << (7 - i));
        }
    }
    return reversed;
}

void NTT(int f1[],long long f_NTT[])
{
        for(int i=0;i<256;i++)
        {
                f_NTT[i]=f1[i];
        }

        long long int start,len,zeta,k=0,tau = 1753,q=8380417;
        long long int t;

        mpz_t q_t, tau_t, zeta_t,t_t,f_t;

        mpz_init(q_t);
        mpz_init(tau_t);
        mpz_init(zeta_t);
        mpz_init(t_t);
        mpz_init(f_t);
        //printf("\n");
        mpz_set_si(tau_t, tau);
        mpz_set_si(q_t, q);
	len = 128;
        //printf("\nZetas in NTT :\n");
	//printf("\n\n values of k: ");

        while(len>=1)
        {
                start = 0;
                while(start < 256)
                {
                        
                        k = k+1;
                        //printf(" k=%d ",k);
                        int r = bitreversal(k);
                        //printf("  reverse=%d  ",r);
                        mpz_powm_ui(zeta_t,tau_t,r,q_t);
                        zeta = mpz_get_si(zeta_t);
                        //zeta=zeta%q;
                        //printf("  zeta=%d  ",zeta);

                        for (int j = start; j<start + len ;j++)
                        {
                                //printf("   zeta=%d  f1[j+len]=%d ",zeta,f1[j+len]);
                                //mpz_set_si(f_t,f1[j+len]);
                                //mpz_mul_si(t_t,f_t,zeta);
                                //t=mpz_get_si(t_t);

                                t = (zeta * f_NTT[j+len])%q;
                                //t=t%q;
                                //printf(" t=%d   ",t);



                                f_NTT[j+len] = (f_NTT[j]-t)%q;
                                if(f_NTT[j+len] < 0)
                                {
                                        f_NTT[j+len]+=q;
                                }

                                f_NTT[j]=(f_NTT[j]+t)%q;
                                if(f_NTT[j]<0)
                                {
                                        f_NTT[j]+=q; 
                                }

                        }
                        start = start + 2* len ;
                }
                len = len/2;
        }

        for(int x = 0; x < 256; x++)
	{
    		if(f_NTT[x] < 0)
    		{
        		f_NTT[x] += q;
    		}
	}

        //printf("\n");

        // printf("\n\nNTT representation of polynomial \n");
        // for(int x = 0; x < 256; x++)
        // {
        //         printf("%d, ",f_NTT[x]);
        // }
        // printf("\n");

}


void INTT(long long int f2[],int q,int tau)
{
        long long int z,start,len,t,zeta;
        mpz_t q_t, tau_t, zeta_t;

        mpz_init(q_t);
        mpz_init(tau_t);
        mpz_init(zeta_t);

        mpz_set_si(tau_t, tau);
        mpz_set_si(q_t, q);
	int k = 256;
	len =1 ;
        //printf("\nZetas in INTT :\n");
        //printf("\n\n values of k: ");
	while (len < 256)
	{
                start = 0;
        	while(start <256)
        	{
                		k = k-1;
                                //printf(" k=%d ",k);

                 		int r = bitreversal(k);
                                //printf("  reverse=%d  ",r);
                        	mpz_powm_ui(zeta_t,tau_t,r,q_t);
                                zeta = mpz_get_si(zeta_t);
                                //printf("  zeta=%d ",zeta);
                                zeta=-zeta%q;
                                if(zeta<0)
                                {
                                        zeta=zeta+q;
                                }
                                //printf(" zeta=%d ",zeta);

                                for ( int j = start ; j <start + len ; j++)
                                {
                                        t = f2[j];
                                        f2[j]=(t + f2[j+len])%q;
					if(f2[j]<0)
                                	{
                                        	f2[j]+=q;
                                	}
                                        f2[j+len]=(t - f2[j+len])%q;
					if(f2[j+len]<0)
                                	{
                                        	f2[j+len]+=q;
                                	}
                                        f2[j+len] = (zeta * f2[j+len])%q;
					if(f2[j+len]<0)
                                	{
                                        	f2[j+len]+=q;
                                	}
                                }
                                start = start + 2*len;
        	}
        	len = 2* len ;

	}
        //printf("\n");

        for(int j = 0; j < 256; j++)
        {
                f2[j]=(f2[j]*8347681)%q;
        }
        
        // for (int i = 0; i < 256; i++)
        // {
        //         out[i]=f2[i];
        // }
        
        // printf("\nINTT representation of polynomial \n");
        // for(int i = 0;i<256;i++)
        // {
        //         printf("%d ",f2[i]);
        // }
        // printf("\n");
}

void multiply_NTT(int a[],long long b[],long long out[])
{
	for (int i = 0; i < 256; i++)
	{
		out[i]=a[i] * b[i];
	}
	
}


void shake_128(unsigned char rho[],unsigned char md_value[],unsigned int md_len,size_t size)
{
		EVP_MD_CTX *mdctx;
     	const EVP_MD *md;
     	    	
     	//unsigned char md_value[EVP_MAX_MD_SIZE];
     	
     	//unsigned int md_len=size;

		int a=size/sizeof(rho[0]);

       	md = EVP_shake128();
 
     	mdctx = EVP_MD_CTX_new();
     	EVP_DigestInit_ex(mdctx, md, NULL);
     	EVP_DigestUpdate(mdctx, rho, a);
		EVP_DigestFinalXOF(mdctx, md_value, md_len);
     	EVP_MD_CTX_free(mdctx);

        //  printf("Md Len = %d ", md_len);
     	// printf("\nDigest is For shake 128 :\n");
     	// for (int i = 0; i < md_len; i++)
     	// {
        //  	printf("%d ,", md_value[i]);
        // }
     	// printf("\n");
     	
     	/*printf("\nDigest is :");
     	for (int i = 0; i < 684; i++)
     	{
         	printf("%02x ", md_value[i]);
        }

     	printf("\n");*/
}


void shake_256(unsigned char input[],unsigned char md_value[],unsigned int md_len,size_t size)
{
		EVP_MD_CTX *mdctx;
     	const EVP_MD *md;
     	    	
     	//unsigned char md_value[EVP_MAX_MD_SIZE];
     	
     	//unsigned int md_len=size;

        int a=size/sizeof(input[0]);

       	md = EVP_shake256();
 
     	mdctx = EVP_MD_CTX_new();
     	EVP_DigestInit_ex(mdctx, md, NULL);
     	EVP_DigestUpdate(mdctx, input, a);
		EVP_DigestFinalXOF(mdctx, md_value, md_len);
     	EVP_MD_CTX_free(mdctx);

        // if(md_len == 8 || md_len == 256){
        // printf("Md Len = %d ", md_len);
     	// printf("\nDigest is For shake 256 :\n");
     	// for (int i = 0; i < md_len; i++)
     	// {
        //  	printf("%d ,", md_value[i]);
        // }
     	// printf("\n");
    // }
}


//key gen


// int coeff_from_three_bytes(int b0,int b1,int b2)
// {
// 	int z;

// 	if(b2 > 127)
// 	{
// 			b2=b2-128;
// 	}
	
// 	z=pow(2,16)*b2+pow(2,8)*b1+b0;
// 	if(z < q)
// 	{
// 		return z;
// 	}
// 	else
// 	{
// 		return -1;			//blank symbol
// 	}
// }

int coeff_from_three_bytes(unsigned char b0, unsigned char b1, unsigned char b2)
{
    int z;

    if (b2 > 127)
    {
        b2 = b2 - 128;
    }

    z = (1 << 16) * b2 + (1 << 8) * b1 + b0;

    if (z < q)
    {
        return z;
    }
    else
    {
        return -1; // blank symbol
    }
}

void RejNTTPoly(unsigned char rho[],int s,int r,int out[])
{
	/*int j_[8],i_[8];
	int_to_bits(j,j_);
	int_to_bits(i,i_);*/

	rho[32]=s;
	rho[33]=r;

	int j=0,c=0;
	
	unsigned char shake128_out[868];
	shake_128(rho,shake128_out,868,34);


	while(j<256)
	{
		out[j]=coeff_from_three_bytes(shake128_out[c],shake128_out[c+1],shake128_out[c+2]);
		c=c+3;
		if(out[j] != -1 )       //blank symbol
		{	
			j++;
		}
	}

    

}


int coeff_from_half_bytes(int b)
{
	//printf("%d  ",b);

	if(eta==2 && b<15)
	{
		return (2-(b % 5));
	}
	else
	{
		if(eta==4 && b<9)
		{
			return 4-b;
		}
		else
		{
			return -9;			//blank symbol
		}
	}
}

// void RejBoundPoly(unsigned char rho_dash[],int out[])
// {
// 	size_t size=66;

// 	int j=0,c=0;
// 	unsigned char shake256_out[512];
// 	shake_256(rho_dash,shake256_out,512,size);
// 	int z,z0,z1; 

// 	/*printf("\nShake out inside RejboundPoly :\n");
// 	for (int i = 0; i < 512; i++)
// 	{
// 		printf("%d ",shake256_out[i]);
// 	}
// 	printf("\n\n");*/
	

// 	while(j<256)
// 	{
// 		z=shake256_out[c];
// 		z0=coeff_from_half_bytes(z%16);
// 		z1=coeff_from_half_bytes(floor(z/16));

// 		if(z0 != -9)
// 		{
// 			out[j]=z0;
// 			j++;
// 		}
// 		if(z1 != -9 && j<256)
// 		{
// 			out[j]=z1;
// 			j++;
// 		}
// 		c++;
// 	}
//     printf("\n\nC value in RejBoundPoly = %d\n",c);
// }

void RejBoundPoly(unsigned char rho_dash[], int out[]) {
    size_t size = 66;  
    unsigned char shake256_out[512];  

    shake_256(rho_dash, shake256_out, 512, size);

    int j = 0, c = 0;
    
    while (j < 256 && c < 512) {  
        unsigned char z = shake256_out[c]; 
        
        unsigned char first_4_bits = z & 0x0F; 
        unsigned char last_4_bits = (z >> 4) & 0x0F;  

        signed char high_bits = (first_4_bits >> 2) & 0x03;  
        signed char low_bits = first_4_bits & 0x03; 

        int coeff1 = (high_bits - low_bits);  

        high_bits = (last_4_bits >> 2) & 0x03; 
        low_bits = last_4_bits & 0x03; 

        int coeff2 = (high_bits - low_bits);

        if (coeff1 >= -2 && coeff1 <= 2 && j < 256) {
            out[j] = coeff1;
            j++;
        }
        if (coeff2 >= -2 && coeff2 <= 2 && j < 256) {
            out[j] = coeff2;
            j++;
        }

        c++;
    }

    printf("\nC value in RejBoundPoly = %d\n", c);
}

void power2Round(int t[],int t0[],int t1[],int d)
{
	int r,r_plus,r0,r1;
	int x=pow(2,d);
	for(int i=0;i<256;i++)
	{
		r=t[i];
		r_plus=r%q;
		r0=( r_plus % x );

		if(r0 > (x/2))
		{
			r0=r0-x;
		}
		if(r0 < -(x/2))
		{
			r0=r0+x;
		}
		
		r1=(r_plus -r0) / x;
		
		t0[i]=r0;
		t1[i]=r1;

		// t1[i] = (t[i] + (1 << (d-1)) - 1) >> d;

		// t0[i] = t[i] - (t1[i] << d);
	}
}


void bits_to_bytes(uint16_t input[],uint16_t out[],int c)    
{
	int x;
	for(int i=0;i<ceil(c/8);i++)
	{
		out[i]=0;
	}
	for(int i=0;i<c;i++)
	{
		x=i%8;
		out[(int)floor(i/8)]=out[(int)floor(i/8)]+input[i]*pow(2,x);
	}

}


void bytes_to_bits(uint16_t input[],uint16_t out[],size_t size)       //output size will be 8*d and d will be length of input
{
	int d=size/sizeof(input[0]);

	uint16_t input1[d]; 
	for(int i=0;i<d;i++)              //
	{
		input1[i]=input[i];
	}

	for(int i=0;i<d;i++)
	{
		for(int j=0;j<8;j++)
		{
			out[8*i+j]=input1[i] % 2;
			input1[i]=floor(input1[i]/2);
		}
	}
}


void int_to_bits(int input,int alpha,uint16_t out[])
{
	int input1=input;
	for(int i=0;i<alpha;i++)
	{
		out[i]=input1%2;
		input1=floor(input1/2);
	}
}

int bits_to_int(uint16_t input[],size_t size)
{
	int x=0;
	int alpha=size;///sizeof(input[0]);
	for(int i=1;i<alpha+1;i++)
	{
		x=2*x+input[alpha-i];
	}
    //printf("%d ",x);
	return x;
}


void simple_bit_pack(int input[],int b,uint16_t out[])
{
	int len=bitLength(b);
    uint16_t z[256*len];
	for(int i=0;i<256*len;i++)
	{
		z[i]=0;
	}

	uint16_t out1[len];
	int count=0;
	for(int i=0;i<256;i++)
	{
		int_to_bits(input[i],len,out1);
		for(int j=0;j<len;j++)
		{
			z[count]=out1[j];
			count++;
		}
	}

	bits_to_bytes(z,out,256*len);
}


void pk_Encode(uint16_t rho[],crystal_poly t1[],uint16_t pk[])
{
	int count=0;
	for(int i=0;i<32;i++)
	{
		pk[i]=rho[i];
		count++;
	}

	int len1=bitLength(q-1);

	uint16_t out[320];
	for(int i=0;i<k;i++)
	{		
		simple_bit_pack( t1[i].arr, ( pow(2,len1-d ) - 1 ) ,out);
		for(int j=0;j<320;j++)
		{
			pk[count]=out[j];
			count++;	
		}
	}
}


void bit_pack(int w[],int a,int b,uint16_t out[])
{
	uint16_t len=bitLength(a+b);
	uint16_t z[len*256],bits[len];

	for (int i = 0; i < len*256; i++)
	{
		z[i]=0;
	}
	

	int count=0;
	for(int i=0;i<256;i++)
	{
		int_to_bits(b-w[i],len,bits);
		for(int j=0;j<len;j++)
		{
			z[count]=bits[j];
			count++;
		}
	}

	bits_to_bytes(z,out,len*256);

}


void sk_encode(uint16_t rho[],unsigned char k1[],unsigned char tr[],crystal_poly s1[],crystal_poly s2[],crystal_poly t0[],uint16_t sk[])
{
	
	int count=0;
	for(int i=0;i<32;i++)
	{
		sk[i]=rho[i];
		count++;
	}
	for(int i=0;i<32;i++)
	{
		sk[count]=k1[i];
		count++;
	}
	for(int i=0;i<64;i++)
	{
		sk[count]=tr[i];
		count++;
	}

	int len=bitLength(eta+eta);
	uint16_t out[len*32];
	for(int i=0;i<l;i++)
	{
		bit_pack(s1[i].arr,eta,eta,out);
		for(int j=0;j<len*32;j++)
		{
			sk[count]=out[j];
			count++;
		}
	}
	for(int i=0;i<k;i++)
	{
		bit_pack(s2[i].arr,eta,eta,out);
		for(int j=0;j<len*32 ;j++)
		{
			sk[count]=out[j];
			count++;
		}
	}

	len=bitLength((pow(2,d-1)-1)+pow(2,d-1));
	uint16_t out1[len*32];
	for(int i=0;i<k;i++)
	{
		bit_pack(t0[i].arr,(pow(2,d-1)-1),pow(2,d-1),out1);
		for(int j=0;j<len*32;j++)
		{
			sk[count]=out1[j];
			count++;
		}
		//printf("\n");
	}
}


//signature

void bit_unpack(uint16_t v[],int out[],int a,int b,size_t size)
{ 
    int c=bitLength(a+b);
    // printf("\nc=%d\n",c);
    // printf("\ny in bit unpack\n");
    // for (int i = 0; i < 32*c; i++)
    // {
    //     printf("%d ",v[i]);
    // }
    // printf("\n");


    uint16_t z[256*c];         //32*c*8
    bytes_to_bits(v,z,size);

    // printf("\nz\n");
    // for (int i = 0; i < 256*c; i++)
    // {
    //     printf("%d ",z[i]);
    // }
    // printf("\n");

    int count;

    uint16_t bits[c];
    
    //printf("\nbits to int op\n");
    for(int i=0;i<256;i++)
    {
        count=0;        
        for(int j=i*c; j<i*c+c ;j++)
        {
            bits[count]=z[j];
            count++;
        }
        out[i]=b - bits_to_int(bits,c);
    }
}


void skDecode(uint16_t sk[],uint16_t rho[],uint16_t k1[],uint16_t tr[],crystal_poly s1[],crystal_poly s2[],crystal_poly t0[],uint16_t eta,uint16_t k,uint16_t l)
{
    // printf("\nsk =\n");
    // for(int i=0;i<4896;i++)
    // {
    //     printf("%d ",sk[i]);
    // }
    // printf("\n");
    uint16_t f[32],g[32],h[64];

    int len=bitLength(2*eta);

    typedef struct arr1
	{
		uint16_t arr[32*len];
	}crystal_poly3;
		
	typedef struct arr2
	{
		uint16_t arr[32*d];
	}crystal_poly4;

    crystal_poly3 y[l],z[k];
	crystal_poly4 w[k];

    int cnt=0;
    for(int i=0;i<32;i++)
    {
        f[i]=sk[i];
        rho[i]=sk[i];
    }
    /*printf("\nRho :");
    for(int i=0;i<32;i++)
    {
        printf("%d ",rho[i]);
    }*/

    for(int i=32;i<64;i++)
    {
        g[cnt]=sk[i];
        k1[cnt]=sk[i];
        cnt++;
    }

    /*printf("\nk1 :");
    for(int i=0;i<32;i++)
    {
        printf("%d ",k1[i]);
    }
    printf("\n");*/

    cnt=0;
    for(int i=64;i<128;i++)
    {
        h[cnt]=sk[i];
        tr[cnt]=sk[i];
        cnt++;
    }


    /*printf("\ntr :");
    for(int i=0;i<64;i++)
    {
        printf("%d ",tr[i]);
    }
    printf("\n");*/

    int size_y=32*len*l,size_z=32*len*k,size_w=32*d*k;
    uint16_t y1[size_y],z1[size_z],w1[size_w];

    cnt=0;
    for(int i=128; i<128+size_y ;i++)
	{
		y1[cnt]=sk[i];
		cnt++;
	}

    cnt=0;
    for(int i=128+size_y; i<128+size_y+size_z ;i++)
	{
		z1[cnt]=sk[i];
		cnt++;
	}

    cnt=0;
    for(int i=128+size_y+size_z; i<128+size_y+size_z+size_w ;i++)
	{
		w1[cnt]=sk[i];
		cnt++;
	}

    cnt=0;
    int count=0;
    for(int i=0;i<l;i++)
    {
        for(int j=i*32*len;count<32*len;j++)
        {
            y[i].arr[cnt]=y1[j];
            cnt++;
            count++;
        }
        cnt=0;
        count=0;
    }

    // printf("\ny :\n");
    // for(int i=0;i<l;i++)
    // {
    //     for(int j=0;j<32*len;j++)
    //     {
    //         printf("%d ",y[i].arr[j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");



    cnt=0;
    count=0;
    for(int i=0;i<k;i++)
    {
        for(int j=i*32*len;count<32*len;j++)
        {
            z[i].arr[cnt]=z1[j];
            cnt++;
            count++;
        }
        cnt=0;
        count=0;
    }

    // printf("\nz :\n");
    // for(int i=0;i<k;i++)
    // {
    //     for(int j=0;j<32*len;j++)
    //     {
    //         printf("%d ",z[i].arr[j]);
    //     }
    // }
    // printf("\n");

    cnt=0;
    count=0;
    for(int i=0;i<k;i++)
    {
        for(int j=i*32*d;count<32*d;j++)
        {
            w[i].arr[cnt]=w1[j];
            cnt++;
            count++;
        }
        cnt=0;
        count=0;
    }

    // printf("\nw :\n");
    // for(int i=0;i<k;i++)
    // {
    //     for(int j=0;j<32*d;j++)
    //     {
    //         printf("%d ",z[i].arr[j]);
    //     }
    // }
    // printf("\n");

    size_t size;
    for(int i=0;i<l;i++)
    {
        size=sizeof(y[i].arr);
        bit_unpack(y[i].arr,s1[i].arr,eta,eta,size);
    }

    // printf("\ns1 inside sk decode:\n");
    // for(int i=0;i<l;i++)
    // {
    //     for(int j=0;j<256;j++)
    //     {
    //         printf("%d ",s1[i].arr[j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    for(int i=0;i<k;i++)
    {
        size=sizeof(z[i].arr);
        bit_unpack(z[i].arr,s2[i].arr,eta,eta,size);
    }

    for(int i=0;i<k;i++)
    {
        size=sizeof(w[i].arr);
        bit_unpack(w[i].arr,t0[i].arr,pow(2,d-1)-1,pow(2,d-1),size);
    }

    //return 0;
}


void ExpandMask(unsigned char rho_dash[],uint16_t neu,crystal_poly y[])
{
    int c=1+bitLength(gamma1 - 1);        
                      
    unsigned char shake256_ip[66],shake256_op[32*c];
    size_t size1=sizeof(shake256_ip);
    int cnt=0;
    int count=0;

    for (int i = 0; i < 64; i++)
    {
        shake256_ip[i]=rho_dash[i];
    }
    
    uint16_t n[16];
    int md_len=32*c;
    uint16_t v[32*c];
    size_t size2=sizeof(v);

    for(int r=0;r<l;r++)
    {
        shake256_ip[64]=neu+r;
        shake256_ip[65]=0;      

        shake_256(shake256_ip,shake256_op,md_len,size1);
        
        // Me- Shake_256 return type is unsigned char, bit_unpack requires uint_16
        count=0;
        for (int i = 0; count < 32*c; i++)
        {
            v[count]=shake256_op[i];
            count++;
        }

        bit_unpack(v,y[r].arr,gamma1 - 1,gamma1,size2); //Me- takes uint_16 and returns int.

    }

}


void decompose(int r,int r2[])
{
    int r1,r0;
    int r_plus=r%q;
    if (r_plus<0)
    {
        r_plus+=q;
    }
    
    
    r0=r_plus % (2*gamma2);
    
    if( r0>gamma2 )
    {
        r0=r0-(2*gamma2);
    }
    else if(r0 < -gamma2)
    {
	    r0=r0+(2*gamma2);
    }

    if((r_plus - r0) == (q-1))
    {
        
        r1=0;
        r0=r0-1;
    }
    else
    {
        r1=(r_plus - r0);
        r1=r1/(2*gamma2);
    }
    
    r2[0]=r1;
    r2[1]=r0;
}


void High_bits(int input[],int out[])
{
    int r[2];
    for (int i = 0; i < 256; i++)
    {
        decompose(input[i],r);
        out[i]=r[0];
    }
    
}

void Low_bits(int input[],int out[])
{
    int r[2];
    for (int i = 0; i < 256; i++)
    {
        decompose(input[i],r);
        out[i]=r[1];
    }

}


void w1_Encode(crystal_poly w1[],uint16_t w1_Encode_out[])
{
    int x=(q-1)/(2*gamma2) - 1;
    int y=bitLength(x);
    int count=0;

    uint16_t out[32*y];
    size_t size=sizeof(out);
    uint16_t bit[256*y];


    for (int i = 0; i < k; i++)
    {
        
        simple_bit_pack(w1[i].arr,x,out);
        
       // bytes_to_bits(out,bit,size);

        for (int j = 0; j < 32*y; j++)
        {
            w1_Encode_out[count]=out[j];
            count++;
        }
        
    }
    
}


// void sample_in_ball(unsigned char c1_bar[],int c_out[],size_t size)
// {
//     // printf("size = %zu\n", size); -> 64
//     printf("c1_bar as hexadecimal:\n");
//     for (size_t i = 0; i < size; i++)
//     {
//         printf("%02x ", c1_bar[i]);  // Print each byte in 2-digit hex format
//     }
//     printf("\n");
    
//     printf("\n Fips Algo cdac \n");
//     int j,c=0,k=0,x;

//     for (int i = 0; i < 256; i++)
//     {
//         c_out[i]=0;
//     }

//     unsigned char shake_256op[64];
//     shake_256(c1_bar,shake_256op,64,size);

//     uint16_t byteArray[8],bitArray[64];
//     for (int i = 0; i < 8; i++)
//     {
//         byteArray[i]=shake_256op[i];
//     }

//     size_t size1=sizeof(byteArray);
//     bytes_to_bits(byteArray,bitArray,size1);

//     for (int i = 256-tau; i < 256; i++)
//     {
//         printf("k = %d, ",k);
//         printf("Shakeop[%d] = %d, ",k,shake_256op[k]);
//         while (shake_256op[k] > i)
//         {
//             printf("%d value exceeded %d ,",shake_256op[k], i);
//             k++;
//         }

//         j=shake_256op[k];
//         c_out[i]=c_out[j];

//         x=i+tau-256;    
//         c_out[j]= pow(-1,bitArray[x]);
//         k++;

//     }

//     printf("\n cout Valyes\n");
//     for(int i=0; i<256; i++){
//         printf("%d ,",c_out[i]);
//     }
//     printf("\n");
    
// }

#define N 256  // Polynomial size
#define HW 60 
void sample_in_ball(unsigned char c1_bar[], int c_out[], size_t size) {
    // printf("c1_bar as hexadecimal:\n");
    // for (size_t i = 0; i < size; i++) {
    //     printf("%02x ", c1_bar[i]);
    // }
    // printf("\n");

    printf("\n HWT Sampling (Hamming Weight = %d)\n", HW);

    int pos = 0;
    uint32_t index, quotient, remainder;
  
    for (int i = 0; i < N; i++) {
        c_out[i] = 0;
    }

    unsigned char shake_256op[256];
    shake_256(c1_bar, shake_256op, 256, size);

    for (int i = 0; i < HW; i++) {

        index = (shake_256op[4 * i + 3] << 24) | (shake_256op[4 * i + 2] << 16) | 
                (shake_256op[4 * i + 1] << 8)  | shake_256op[4 * i];

        quotient = 0xffffffff / (N - HW + pos + 1);
        remainder = 0xffffffff - (quotient * (N - HW + pos + 1));

        // printf("Index = %u || ", index);
        // printf("Remainder = %u || ", remainder);
        // printf("Quotient = %u || ", quotient);

        if ((0xffffffff - index) > remainder && pos < HW) {
            index = index / quotient;
            printf("Index/q = %u || ", index);

            if (index >= N) continue; 

            c_out[N - HW + pos] = c_out[index];

            int bit_pos = HW * 4 + (i >> 3);
            if (bit_pos < sizeof(shake_256op)) {
                c_out[index] = ((((shake_256op[bit_pos] >> (i & 0x07)) & 0x01) << 1) - 1);
            }

            pos++;
        }
    }
    
    if (pos != HW) {
        printf("Error: Hamming weight mismatch! Expected %d, got %d\n", HW, pos);
    }

    printf("\nc_out Values:\n");
    for (int i = 0; i < N; i++) {
        printf("%d, ", c_out[i]);
    }
    printf("\n");
}

 
// void sample_in_ball(unsigned char c1_bar[], int c_out[], size_t size) {

//     int indices[N];
//     unsigned char shake_256op[N];

//     for (int i = 0; i < N; i++) {
//         c_out[i] = 0;
//         indices[i] = i;
//     }

//     shake_256(c1_bar, shake_256op, N, size);

//     printf("\nGenerated shake_256op values:\n");
//     for (int i = 0; i < 256; i++) {
//         printf("[%d] -> %d ,", i,shake_256op[i]);
//     }
//     printf("\n");

//     for (int i = N - 1; i >= N - TAU; i--) {
//         int swap_idx = shake_256op[i] % (i + 1); 
//         int temp = indices[i];
//         indices[i] = indices[swap_idx];
//         indices[swap_idx] = temp;
//     }

//     printf("\n Indices \n");
//     for(int i=0; i<256; i++){
//         printf("[%d] -> %d ,",i,indices[i]);
//     }
//     printf("\n");

//     for (int i = N - TAU; i < N; i++) {
//         int pos = indices[i];
//         c_out[pos] = (shake_256op[i] & 1) ? 1 : -1; 
//     }

//     printf("\n C out: \n");
//     for(int i=0; i<256; i++){
//         printf("[%d] -> %d ,",i,c_out[i]);
//     }
//     printf("\n");
// }

int inf_norm_Grp(int z)
{
//int q = 17;
z = z % q;


if(z>(q/2))
{
	z=z-q;
}
if(z<(-q/2))
{
	z= z+q;
}
return z;
}


int inf_norm_poly(int z[])
{
	int out[256];
	for(int i =0;i<256;i++)
	{
		out[i]=inf_norm_Grp(z[i]);
        if(out[i]<0)
        {
            out[i]=abs(out[i]);
        }
	}

	int max = out[0] ;
    for (int i=1;i<256;i++) 
    { 
        if ( out[i] > max )
        {
            max = out[i];
        }             
    }

    return max ;  	 
}


int inf_norm_poly_ct0(long long int z[])
{
	int out[256];
	for(int i =0;i<256;i++)
	{
		out[i]=inf_norm_Grp(z[i]);
        if(out[i]<0)
        {
            out[i]=abs(out[i]);
        }
	}

	int max = out[0] ;
    for (int i=1;i<256;i++) 
    { 
        if ( out[i] > max )
        {  
            max = out[i] ;  
        }
    }
    
    return max ;  	 
}


int make_hint(crystal_poly1 z[],crystal_poly r[],crystal_poly h[])
{
    crystal_poly r1[k],v1[k],input[k];
    int count=0;

    for (int i = 0; i < k; i++)
    {
        High_bits(r[i].arr,r1[i].arr);
    }
    
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < 256; j++)
        {
            input[i].arr[j]=r[i].arr[j] + z[i].arr[j];
        }
    }

    // for (int i = 0; i < k; i++)
    // {
    //     High_bits(r[i].arr,r1[i].arr);
    // }

    for (int i = 0; i < k; i++)
    {
        High_bits(input[i].arr,v1[i].arr);
    }

    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < 256; j++)
        {
            if(r1[i].arr[j] != v1[i].arr[j])
            {
                h[i].arr[j] = 1;
                //count++;
            }
            else
            {
                h[i].arr[j] = 0;
            }
        }
    }

    

    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < 256; j++)
        {
            if(h[i].arr[j] == 1)
            {             
                count++;
            }
        }
    }

    return count;

}


void hint_bit_pack(crystal_poly h[],int y[])
{
    int index=0;

    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < 256; j++)
        {
            if(h[i].arr[j] != 0)
            {
                y[index]=j;
                index++;
            }
        }
        y[omega+i]=index;
    } 
}


void sig_encode(unsigned char c_bar[],crystal_poly z[],crystal_poly h[],uint16_t sign[])
{
    int count=0;
    for (int i = 0; i < (2*lamda)/8; i++)
    {
        sign[i]=c_bar[i];
        count++;
    }
    
    int len=bitLength(gamma1-1+gamma1);
    uint16_t out[32*len];
    for (int i = 0; i < l; i++)
    {
        bit_pack(z[i].arr,gamma1-1,gamma1,out);
        for (int j = 0; j < 32*len; j++)
        {
            sign[count]=out[j];
            count++;
        }
    }

    int y[omega+k];
    for (int i = 0; i < omega+k; i++)
    {
        y[i]=0;
    }
    

    hint_bit_pack(h,y);
    for (int i = 0; i < omega+k; i++)
    {
        sign[count]=y[i];
        count++;
    }
        
}


//verification

void simple_bit_unpack(uint16_t v[],int out[],int b,size_t size)
{
    int c=bitLength(b);

    uint16_t z[256*c];         //32*c*8
    bytes_to_bits(v,z,size);

    int count;

    uint16_t bits[c];
    
    for(int i=0;i<256;i++)
    {
        count=0;        
        for(int j=i*c; j<i*c+c ;j++)
        {
            bits[count]=z[j];
            count++;
        }
        out[i]=bits_to_int(bits,c);
    }


}


void pk_decode(uint16_t pk[],uint16_t rho[],crystal_poly t1[])
{
    for (int i = 0; i < 32; i++)
    {
        rho[i]=pk[i];
    }

    int len=32*(bitLength(q-1)-d);

    typedef struct arr3
	{
		uint16_t arr[len];
	}crystal_poly5;

    uint16_t z1[k*len];

    crystal_poly5 z[k];

    int count=0;
    for (int i = 32; i < 32+k*len; i++)
    {
        z1[count]=pk[i];
        count++;
    }

    int cnt=0;
    for (int i = 0; i < k; i++)
    {
        cnt=0;
        for (int j = i*len; cnt < len; j++)
        {
            z[i].arr[cnt]=z1[j];
            cnt++;
        }    
    }
    
    size_t size;
    for(int i=0;i<k;i++)
    {
        size=sizeof(z[i].arr);
        simple_bit_unpack(z[i].arr,t1[i].arr,pow(2,(bitLength(q-1)-d)) -1 ,size);
    }

}


int hint_bit_unpack(uint16_t y[],crystal_poly h[])
{
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < 256; j++)
        {
            h[i].arr[j]=0;
        }
    }

    int index=0,x;
    for (int i = 0; i < k; i++)
    {
        if(y[omega+i]<index || y[omega+i]>omega)
        {
            return -1;                      //rejection symbol
        }

        while(index<y[omega+i])
        {
            x=y[index];
            h[i].arr[x]=1;
            index++;
        }
    }

    while (index<omega)
    {
        if(y[index] != 0)
        {
            return -1;                     //rejection symbol
        }
        index++;
    }
   
}


// void sig_decode(uint16_t sign[],unsigned char c_bar[],crystal_poly z[],crystal_poly h[])
// {
    
//     for (int i = 0; i < 32; i++)
//     {
//         c_bar[i]=sign[i];
//     }

//     int len1=32*(1+bitLength(gamma1-1));
//     uint16_t x1[l*len1];

//     typedef struct arr4
// 	{
// 		uint16_t arr[len1];
// 	}crystal_poly6;

//     crystal_poly6 x[l];

//     int count=0;
//     for (int i = 32; i < 32 + len1*l; i++)
//     {
//         x1[count]=sign[i];
//         count++;
//     }
    
//     // printf("\nx1 =\n");
//     // for (int i = 0; i < l*len1; i++)
//     // {
//     //     printf("%d ",x1[i]);
//     // }
//     // printf("\n");

//     count=0;
//     for (int i = 0; i < l; i++)
//     {
//         count=0;
//         for (int j = i*len1; count < len1; j++)
//         {
//             x[i].arr[count]=x1[j];
//             count++;
//         }    
//     }

//     // printf("\nx =\n");
//     // for (int i = 0; i < l; i++)
//     // {
//     //     for (int j = 0; j < len1; j++)
//     //     {
//     //         printf("%d ",x[i].arr[j]);
//     //     }
//     // }
//     // printf("\n");
    
    
//     size_t size;
//     for(int i=0;i<l;i++)
//     {
//         size=sizeof(x[i].arr);
//         bit_unpack(x[i].arr,z[i].arr,gamma1-1,gamma1,size);
//     }

    
//     uint16_t y[omega+k];
//     count=0;
//     for (int i = (l*len1)+32; i < 32+(l*len1)+omega+k; i++)
//     {
//         y[count]=sign[i];
//         count++;
//     }
    
//     int abc=hint_bit_unpack(y,h);
//     //printf("\nafter bit hint bit unpack\n");
// }




void sig_decode(uint16_t sign[],unsigned char c_bar[],crystal_poly z[],crystal_poly h[])
{
    for (int i = 0; i < (2*lamda)/8; i++)
    {
        c_bar[i]=sign[i];
    }

    int len1=32*(1+bitLength(gamma1-1));
    uint16_t x1[l*len1];

    typedef struct arr4
	{
		uint16_t arr[len1];
	}crystal_poly6;

    crystal_poly6 x[l];

    int count=0;
    for (int i = (2*lamda)/8; i < ((2*lamda)/8)+len1*l; i++)
    {
        x1[count]=sign[i];
        count++;
    }
    
    count=0;
    for (int i = 0; i < l; i++)
    {
        count=0;
        for (int j = i*len1; count < len1; j++)
        {
            x[i].arr[count]=x1[j];
            count++;
        }    
    }

    size_t size;
    for(int i=0;i<l;i++)
    {
        size=sizeof(x[i].arr);
        bit_unpack(x[i].arr,z[i].arr,gamma1-1,gamma1,size);
    }

    uint16_t y[omega+k];
    count=0;
    for (int i = l*len1+((2*lamda)/8); i < ((2*lamda)/8)+l*len1+omega+k; i++)
    {
        y[count]=sign[i];
        count++;
    }
    
    hint_bit_unpack(y,h);
    
}


// void use_hint(int h[],long long r[],int w1[])
// {
//     int m=(q-1)/(2*gamma2);
//     int r2[2];
   
//     for (int i = 0; i < 256; i++)
//     {
//         decompose(r[i],r2);

//         //printf("%d ",r2[0]);

//         if(h[i]==1 && r2[1]>0)
//         {
//             w1[i]=(r2[0]+1) %m;//& 15;
//         }
//         if(h[i]==1 && r2[1]<=0)
//         {
            
//             if(r2[0]==0 && gamma2==95232)
//             {
//                 w1[i]=43;
//             }
//             else
//             {
//                 w1[i]=(r2[0]-1) & 15;
//             }
//         }
//         if(h[i]==0)
//         {
//             w1[i]=r2[0];
//         }
//     }
// }



void use_hint(int h[],long long r[],int w1[])
{
    int m=(q-1)/(2*gamma2);
    int r2[2];
   
    for (int i = 0; i < 256; i++)
    {
        decompose(r[i],r2);
        if(h[i]==1 && r2[1]>0)
        {
            w1[i]=(r2[0]+1)%m;
        }
        if(h[i]==1 && r2[1]<=0)
        {
            w1[i]=(r2[0]-1)%m;
            if(w1[i]<0)
            {
                w1[i]=w1[i]+m;
            }
        }
        if(h[i]==0)
        {
            w1[i]=r2[0];
        }
    }
}