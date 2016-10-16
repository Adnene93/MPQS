//compile with
//gcc -o factorize factorize.c -lgmp -lm

//mpirun -np 1 --hostfile=hostfile ./factorize n b m 

#include <stdio.h>
#include <stdlib.h>			/* for printf */
#include <math.h>
#include <gmp.h>
#include <time.h>
#include <omp.h>
#include "mpi.h"
#define MAXVECTOR_SIZE 100000
#define MAXVECTORHORSBASE_SIZE 1000
#define MAXPOLYNOMES 50000

//#define verbose 1

#define get_elem2(indexprime,index,nbit) \
	index = indexprime/8;\
	nbit = indexprime%8;


unsigned int verbose = 1;

int nbprocess,rank;
int *stops;
//standard
int prime_test(mpz_t number,long iter)
{
    int returned = 1;
    if ((mpz_cmp_ui(number,2)==0) || (mpz_cmp_ui(number,3)==0)) return 1;
    mpz_t decremented;
    mpz_init(decremented);
    mpz_sub_ui(decremented,number,1); //decremented <- n-1
    mpz_t d;
    mpz_init(d);
    mpz_set(d,decremented); //d <- decremented
    long s = mpz_scan1(decremented,0); //
    if (s==0) return 0; //number is pair
    long i=0;
    int gotothenextloop;
    while (i<s)
    {
        mpz_cdiv_q_ui(d,d,2);
        i++;
    }
    //decremented = 2^s * d
    mpz_t a;
    mpz_init(a);
    mpz_t x;
    mpz_init(x);
    mpz_t max_random;
    mpz_init(max_random);
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    mpz_sub_ui(max_random,number,3); //max = n-3
    while (iter>0)
    {
        mpz_urandomm (a, state,max_random); // a is random and 2 <= a <= n-2
        mpz_add_ui(a,a,2);
        mpz_powm(x,a,d,number); //x = a^d mod n
        if (((mpz_cmp_ui(x,1)==0) || (mpz_cmp(x,decremented)==0))==0)
        {
            i=0;
            gotothenextloop=0;
            while (i<s-1 && !gotothenextloop)
            {
                mpz_powm_ui (x, x, 2, number);
                if (mpz_cmp_ui(x,1)==0) return 0;
                if (mpz_cmp(x,decremented)==0) gotothenextloop=1;
                i++;
            }
            if (!gotothenextloop) return 0;
        }
        iter--;
    }
    return 1;
}

//upgraded
int prime_test2(mpz_t number,long iter)
{
    int returned = 1;
    unsigned long firstprimes[170] = {2,3,5,7,11,13,17,19,23,29, 31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,
                                      127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,    227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,    331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,    439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,    563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,    673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,
                                      811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,
                                      937,941,947,953,967,971,977,983,991,997,1009,1013
                                     };


    unsigned long arrayofPSWJ[12]=     {2047,1373653,9080191,25326001,4759123141,
                                        1122004669633,2152302898747,3474749660383,341550071728321,
                                        3825123056546413051
                                       };
    unsigned long arrayoftesters[12][10] = {{1,2},{2,2,3},{2,31,37},{3,2,3,5},{3,2,7,61},{4,2,13,23,1662803},{5,2,3,5,7,11},{6,2,3,5,7,11,13},{7,2,3,5,7,11,13,17},{9,2, 3, 5, 7, 11, 13, 17, 19,23}};
    long i=0,j=0;

    if (mpz_cmp_ui(number,1014)<0)
    {
        for (i=0; i<170; i++)
        {
            if (mpz_cmp_ui(number,firstprimes[i])==0) return 1;
        }
        return 0;
    }

    mpz_t decremented;
    mpz_init(decremented);
    mpz_sub_ui(decremented,number,1);
    mpz_t d;
    mpz_init(d);
    mpz_set(d,decremented);
    long s = mpz_scan1(decremented,0);
    if (s==0) return 0;
    int gotothenextloop;
    i=0;
    while (i<s)
    {
        mpz_cdiv_q_ui(d,d,2);
        i++;
    }

    mpz_t a;
    mpz_init(a);
    mpz_t x;
    mpz_init(x);
    mpz_t max_random;
    mpz_init(max_random);
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    mpz_sub_ui(max_random,number,3);
    i=0;
    int k=0;
    for (k=0; k<10; k++)
    {
        if (mpz_cmp_ui(number,arrayofPSWJ[k])<0) break;
    }

    if (k<10)
    {
        i=1;
        for (j=1; j<=arrayoftesters[k][0]; j++)
        {
            mpz_set_ui(a,arrayoftesters[k][j]);
            mpz_powm(x,a,d,number);
            if (((mpz_cmp_ui(x,1)==0) || (mpz_cmp(x,decremented)==0))==0)
            {
                i=0;
                gotothenextloop=0;
                while (i<s-1 && !gotothenextloop)
                {
                    mpz_powm_ui (x, x, 2, number);
                    if (mpz_cmp_ui(x,1)==0) return 0;
                    if (mpz_cmp(x,decremented)==0) gotothenextloop=1;
                    i++;
                }
                if (!gotothenextloop) return 0;
            }
        }
    }
    else
    {
        while (iter>0)
        {
            mpz_urandomm (a, state,max_random);
            mpz_add_ui(a,a,2);
            mpz_powm(x,a,d,number);
            if (((mpz_cmp_ui(x,1)==0) || (mpz_cmp(x,decremented)==0))==0)
            {
                i=0;
                gotothenextloop=0;
                while (i<s-1 && !gotothenextloop)
                {
                    mpz_powm_ui (x, x, 2, number);
                    if (mpz_cmp_ui(x,1)==0) return 0;
                    if (mpz_cmp(x,decremented)==0) gotothenextloop=1;
                    i++;
                }
                if (!gotothenextloop) return 0;
            }
            iter--;
        }
    }
    return 1;
}



void factorize_semiprime(mpz_t number,mpz_t *p,mpz_t *q)
{
    mpz_t rechp,rechq;
    mpz_init(rechp);
    mpz_init(rechq);
    mpz_sqrt(rechp, number);
    mpz_t quotient,remainder;
    mpz_init(quotient);
    mpz_init(remainder);

    mpz_tdiv_qr_ui (quotient,remainder, rechp, 6);
    int decrementer=0;
    if (mpz_cmp_ui(remainder,0)==0)
    {
        mpz_sub_ui(rechp,rechp,1);
        decrementer=4;
    }
    else if (mpz_cmp_ui(remainder,1)==0) decrementer=2;
    else if (mpz_cmp_ui(remainder,2)==0)
    {
        mpz_sub_ui(rechp,rechp,1);
        decrementer=2;
    }
    else if (mpz_cmp_ui(remainder,3)==0)
    {
        mpz_sub_ui(rechp,rechp,2);
        decrementer=2;
    }
    else if (mpz_cmp_ui(remainder,4)==0)
    {
        mpz_sub_ui(rechp,rechp,3);
        decrementer=2;
    }
    else if (mpz_cmp_ui(remainder,5)==0)
    {
        //mpz_sub_ui(rechp,rechp,3);
        decrementer=4;
    }
    /*printf("%s = %s * 6 + %s | %d\n",mpz_get_str (NULL, 10, rechp),mpz_get_str (NULL, 10, quotient),mpz_get_str (NULL, 10, remainder),decrementer);*/

    printf("Number :%s, start : %s\n",mpz_get_str (NULL, 10, number),mpz_get_str (NULL, 10, rechp));
    mpz_set_ui(remainder,1);
    while (mpz_cmp_ui(remainder,0)!=0)
    {
        printf("rechp = %s\n",mpz_get_str (NULL, 10, rechp));
        //if (prime_test2(rechp,20)==1) {
        mpz_tdiv_qr(rechq, remainder,number, rechp);
        //printf("isprime\n");
        //}
        if (mpz_cmp_ui(remainder,0)!=0) mpz_sub_ui(rechp,rechp,decrementer);
        if (decrementer==2) decrementer=4;
        else decrementer=2;
        //getchar();
    }
    /*printf("%s = %s * %s\n",mpz_get_str (NULL, 10, number),mpz_get_str (NULL, 10,     rechp),mpz_get_str (NULL, 10, rechq));*/
    mpz_set(*p,rechp);
    mpz_set(*q,rechq);
}

void factorize_semiprimegoldbach(mpz_t number,mpz_t *p,mpz_t *q)
{
    mpz_t n,m;
    mpz_init(n);
    mpz_init(m);
    mpz_t squaren,squarem;
    mpz_init(squaren);
    mpz_init(squarem);
    mpz_sqrt(n, number);
    mpz_add_ui(n,n,1);
    mpz_pow_ui(squaren,n,2);
    mpz_sub(squarem,squaren,number);
    while (!mpz_perfect_square_p(squarem))
    {
        mpz_add_ui(n,n,1);
        mpz_pow_ui(squaren,n,2);
        mpz_sub(squarem,squaren,number);
    }

    /*mpz_t rechp,rechq;
    mpz_init(rechp);mpz_init(rechq);
    mpz_t quotient,remainder;
    mpz_init(quotient);mpz_init(remainder);*/

    mpz_sqrt(m,squarem);
    mpz_set(*p,n);
    mpz_sub(*p,*p,m);
    mpz_set(*q,n);
    mpz_add(*q,*q,m);


}

void factorize_semiprime2(mpz_t number,mpz_t *p,mpz_t *q)
{
    mpz_t x,d;
    mpz_init(x);
    mpz_init(d);
    if (mpz_even_p (number))
    {
        mpz_set_ui(*p,2);
        mpz_tdiv_q(*q,number, *p);
    }
    else if (mpz_perfect_square_p(number))
    {
        mpz_sqrt(*p,number);
        mpz_set(*q,*p);
    }
    else
    {
        mpz_t squarex,squared;
        mpz_init(squarex);
        mpz_init(squared);
        mpz_sqrt(x, number);
        mpz_add_ui(x,x,1);
        mpz_t divider;
        mpz_init(divider);
        mpz_t divider2;
        mpz_init(divider2);
        mpz_set_ui(divider2,1);

        while ((mpz_cmp_ui(divider2,1)==0) || (mpz_cmp(divider2,number)==0))
        {
            mpz_powm_ui(squarex,x,2,number);
            while (!mpz_perfect_square_p(squarex))
            {
                mpz_add_ui(x,x,1);
                mpz_powm_ui(squarex,x,2,number);
            }
            mpz_sqrt(d,squarex);
            mpz_add(divider,x,d);
            mpz_gcd(divider2,number,divider);
        }

        mpz_set(*p,divider);
        mpz_tdiv_q(*q,number, divider);
    }
    //mpz_set(*q,divider2);
}

inline void next_prime(mpz_t prime,mpz_t *next)
{
    mpz_t r;
    mpz_init(r);
    mpz_t temp;
    mpz_init(temp);
    mpz_set(temp,prime);
    while (!prime_test2(temp,18))
    {
        mpz_sub_ui(temp,temp,1);
    }
    if (mpz_cmp_ui(temp,2)==0) mpz_set_ui(*next,3);
    else if (mpz_cmp_ui(temp,3)==0) mpz_set_ui(*next,5);
    else
    {

        mpz_mod_ui(r,temp,6);
        int inc;

        if (mpz_cmp_ui(r,1)==0) inc=4;
        else inc=2;
        mpz_set(*next,temp);
        mpz_add_ui(*next,*next,inc);
        if (inc==2) inc=4;
        else inc =2;
        while (!prime_test2(*next,18))
        {
            mpz_add_ui(*next,*next,inc);
            if (inc==2) inc=4;
            else inc =2;
        }
    }
    mpz_clear(r);
    mpz_clear(temp);

}
inline int isQuadraticResidue(mpz_t a,mpz_t modulo)  //modulo is prime
{
    mpz_t puissance,apowered;
    mpz_init(puissance);
    mpz_init(apowered);
    mpz_sub_ui(puissance,modulo,1);
    mpz_tdiv_q_ui(puissance,puissance,2);
    mpz_powm(apowered,a,puissance,modulo);
    //printf("%s\n",mpz_get_str (NULL, 10,apowered));
    if (mpz_cmp_ui(apowered,1) == 0 )
    {
        mpz_clear(apowered);
        mpz_clear(puissance);
        return 1;
    }
    mpz_clear(apowered);
    mpz_clear(puissance);
    return 0;
}

//Tonelli–Shanks_algorithm //a must be coprime with modulo
inline int findQuadraticResidue(mpz_t a,mpz_t modulo,mpz_t *residue1,mpz_t *residue2)
{
    if (!isQuadraticResidue(a,modulo))
    {
        mpz_set_ui(*residue1,0);
        mpz_set_ui(*residue2,0);
        return 0;
    }
    else
    {
        if (mpz_cmp_ui(modulo,2)==0)
        {

            if (mpz_even_p (a)) mpz_set_ui(*residue1,0); //solution triviale
            else mpz_set_ui(*residue1,1);
            mpz_set_ui(*residue2,0);
            return 1;
        }
    }
    mpz_t decremented;
    mpz_init(decremented);
    mpz_sub_ui(decremented,modulo,1); //decremented <- p-1 where p = modulo
    mpz_t q;
    mpz_init(q);
    mpz_set(q,decremented); //q <- p-1
    long s = mpz_scan1(decremented,0); //
    long i=0;
    while (i<s)
    {
        mpz_cdiv_q_ui(q,q,2);
        i++;
    }
    //p-1 = 2^s * q

    mpz_t power;
    mpz_init(power);
    mpz_t z;
    mpz_init(z);


    if (s == 1)
    {
        mpz_add_ui(power,modulo,1);
        mpz_tdiv_q_ui(power,power,4);
        mpz_powm(*residue1,a,power,modulo); //r1 = a^((p+1)/4) mod p
        mpz_sub(*residue2,modulo,*residue1);
        mpz_mod(*residue2,*residue2,modulo); //r2 = p - r1 mod p
    }
    else
    {

        mpz_set_ui(z,2);
        while (isQuadraticResidue(z,modulo)) mpz_add_ui(z,z,1);

        mpz_powm(z,z,q,modulo); // z = z^q mod p

        mpz_add_ui(power,q,1);
        mpz_tdiv_q_ui(power,power,2); // power = (Q+1)/2
        mpz_powm(*residue1,a,power,modulo); // R = a^((q+1)/2)

        mpz_t t;
        mpz_init(t);
        mpz_powm(t,a,q,modulo); // t = a^q mod p
        long m=s;
        mpz_t b;
        mpz_init(b);
        i=1;
        long j = 0;
        mpz_t t2;
        mpz_init(t2);
        while (mpz_cmp_ui(t,1)!=0)
        {
            i=1;
            j=0;
            mpz_set(t2,t);
            while (i<m && mpz_cmp_ui(t2,1)!=0)
            {
                mpz_powm_ui(t2,t2,2,modulo);
                if (mpz_cmp_ui(t2,1)==0) break;
                i++;
            }
            mpz_set(b,z); // b = z;
            while (j<m-i-1)
            {
                mpz_powm_ui(b,b,2,modulo);
                j++;
            } // b = Z^(2^(m-i-1))

            mpz_mul(*residue1,*residue1,b);
            mpz_mod(*residue1,*residue1,modulo);

            mpz_powm_ui(z,b,2,modulo);

            mpz_mul(t,t,z);

            mpz_mod(t,t,modulo);
            m = i;

        }

        mpz_sub(*residue2,modulo,*residue1);
        mpz_mod(*residue2,*residue2,modulo); //r2 = p - r1 mod p

    }
    mpz_clear(z);
    mpz_clear(q);
    mpz_clear(power);
    mpz_clear(decremented);
    return 1;
}

inline void findFirstPrimesVerifieQuadraticResidue(mpz_t n,mpz_t array[],unsigned long size)  //which n-first primes where n have quadratic residue
{
    long i=0;
    for (i=0; i<size; i++)
    {
        mpz_init((array)[i]);
    }
    mpz_t prime,next;
    mpz_init(prime);
    mpz_init(next);
    mpz_set_ui(prime,2);
    mpz_set_ui((array)[0],1); //1 Represent -1
    i=1;
    //i=0;
    while (i<size)
    {
        if (isQuadraticResidue(n,prime))
        {

            mpz_set((array)[i],prime);
            i++;
        }
        next_prime(prime,&next);
        mpz_set(prime,next);
    }

}

inline void findFirstPrimesVerifieQuadraticResidue2(mpz_t n,unsigned int array[],unsigned long size)  //which n-first primes where n have quadratic residue
{
    long i=0;
    mpz_t prime,next;
    mpz_init(prime);
    mpz_init(next);
    mpz_set_ui(prime,2);
    //mpz_set_ui((array)[0],1); //1 Represent -1
    array[0]=1;
    i=1;
    //i=0;
    while (i<size)
    {
        if (isQuadraticResidue(n,prime))
        {

            //mpz_set((array)[i],prime);
            array[i] = mpz_get_ui(prime);
            i++;
        }
        mpz_nextprime(prime,prime);
        //next_prime(prime,&next);
        //mpz_set(prime,next);
        //mpz_nextprime (prime, prime);
    }
    mpz_clear(prime);
    mpz_clear(next);

}

inline void findFirstPrimesVerifieQuadraticResidue3(mpz_t n,unsigned int array[],unsigned long bound,unsigned long *size)  //which n-first primes where n have quadratic residue
{
    long i=0;
    *size=0;
    mpz_t prime,next;
    unsigned long primelong=0;
    mpz_init(prime);
    mpz_init(next);
    mpz_set_ui(prime,2);
    //mpz_set_ui((array)[0],1); //1 Represent -1
    array[0]=1;
    i=1;
    //i=0;
    while (primelong<bound)
    {
        if (isQuadraticResidue(n,prime))
        {

            //mpz_set((array)[i],prime);
            array[*size] = mpz_get_ui(prime);
            (*size) = (*size) + 1;
        }
        next_prime(prime,&next);
        mpz_set(prime,next);
        primelong = mpz_get_ui(prime);
    }
    mpz_clear(prime);
    mpz_clear(next);

}


inline int isNullVector(unsigned char *array,int size)
{
    int i=0;
    int or=0;
    for (i=0; i<size; i++)
    {
        if ((i==0) && (array[i]==1))
        {
            or=0;
        }
        else or+=array[i];
        if (or!=0) return 0;
    }
    return (or==0);
}

unsigned char xor(unsigned char a,unsigned char b)
{
    //if (((a==0) && (b==0)) || ((a==1) && (b==1))) return 0;
    //else return 1;
    unsigned char c;


    return a^b;
}


inline void get_elem(unsigned long indexprime,unsigned int *index,unsigned int *nbit)
{
    (*index) = indexprime/8;
    (*nbit) = indexprime%8;
}

unsigned char get_elem_value(unsigned long i,unsigned long indexprime,unsigned char **matrix)
{
    unsigned int index = indexprime/8;
    unsigned int  nbit = indexprime%8;
    unsigned char x =  (matrix[i][index] & (1 << nbit)) >> nbit;
    //printf("\nx=%d\n",x);
    if((matrix[i][index] & (1 << nbit)) >> nbit) return 1;
    else return 0;
}

inline double logarithm10(mpz_t n){
        mpf_t x,tmp,res;
        mpf_init(x);
        mpf_init(tmp);
        mpf_init(res);
        mpf_set_z(x,n);
        int b = (mpz_sizeinbase (n, 10));
        mpf_set_ui(tmp,10);
        mpf_pow_ui(tmp,tmp,b);
        mpf_div(res,x,tmp);
        mpf_sub_ui(res,res,1);
        
        mpf_set(x,res);
        int i=2;
        unsigned char add = 0;
        while (i<100){
            mpf_set(tmp,x);
            mpf_pow_ui(tmp,tmp,i);
            mpf_div_ui(tmp,tmp,i);
            if (!add) {mpf_sub(res,res,tmp);add=1;}
            else {mpf_add(res,res,tmp);add=0;}
            i+=1;
        }
     
        double ret = mpf_get_d(res);
        ret+=((double) b)*2.302585093;
        ret = ret/2.302585093;
        printf("\nlog = %f\n",ret); 
	return ret;
}


/*inline void next_prime_QuadraticResidue(mpz_t n,mpz_t prime,mpz_t *next)

{

    	next_prime(prime,next);
        while (!isQuadraticResidue(n,*next)) next_prime(*next,next);        

}*/

inline void next_prime_QuadraticResidue(mpz_t n,mpz_t prime,mpz_t *next)

{

    	mpz_nextprime(*next,prime);
        while (!isQuadraticResidue(n,*next)) mpz_nextprime(*next,*next);        

}


inline void next_k_prime_QuadraticResidue(mpz_t n,mpz_t prime,mpz_t *next,int k)

{
    int boucle = 0;
    mpz_set(*next,prime);

    for (boucle=0;boucle<k;boucle++) next_prime_QuadraticResidue(n,*next,next);

            

}
void factorize(mpz_t n,mpz_t *p,mpz_t *q,long nprimes,long poolsize)
{
    
    double timebeginseive=0;
    double timeendseive=0;
    unsigned long maxsmooths=nprimes+nprimes/4;
    unsigned int numprime=0,numbit=0,numprime2=0,numbit2=0;
    unsigned long numpoly=0;
    unsigned char shouldclear=0;
    unsigned long i=0,j=0,j2=0,i1=0,i22=0,maxalloc=0;
    if (maxsmooths<1000) maxalloc=1000;
    else maxalloc=maxsmooths+nprimes/4;
    unsigned int counting2;
    unsigned long k;
    char  *s;
    unsigned long inc1,inc2,actualincrementer;
    long test,test2;
    s = malloc((2015)*sizeof(char));
    
    mpz_t rootnsquare;
    mpz_init(rootnsquare);
    mpz_t rootnprime;
    mpz_init(rootnprime);
    mpz_t tmp;
    mpz_init(tmp);
    mpz_t prime;
    mpz_init(prime);
    mpz_t r1,r2;
    mpz_init(r1);
    mpz_init(r2);
    mpz_t rootn;
    mpz_init(rootn);
    mpz_sqrt(rootn, n);
    mpz_mul(rootnsquare,rootn,rootn);
    if (prime_test2(n,18)) {
	if (rank==0) printf("this number is prime\n"); 
        mpz_set(*p,n);
        mpz_set_ui(*q,1);
        return;
    } 
    if (mpz_cmp(rootnsquare,n)==0)
    {
        mpz_init(*p);
        mpz_init(*q);
        mpz_set(*p,rootn);
        mpz_set(*q,rootn);
        if (verbose && rank==0) printf("\nLucky : This Is A Perfect Square");
        return;
    }

    mpz_t A,invA,b,B,Btemp,C,m_val;
    mpz_init(A);
    mpz_init(invA);
    mpz_init(b);
    mpz_init(B);
    mpz_init(Btemp);
    mpz_init(C);
    mpz_init(m_val);
    //----------rootnprime = racine(racine(2*n)/M)---------
    mpz_set(rootnprime,n);
    mpz_mul_ui(rootnprime,rootnprime,2);
    mpz_sqrt(rootnprime,rootnprime);
    mpz_tdiv_q_ui(rootnprime,rootnprime,poolsize);
    mpz_sqrt(rootnprime,rootnprime);
    if (mpz_cmp_ui(rootnprime,1)==0) mpz_set_ui(rootnprime,2);
    mpz_t rootnprimebegin;
    mpz_init(rootnprimebegin);
    mpz_set(rootnprimebegin,rootnprime);
    mpz_add_ui(rootn,rootn,1);
    //----------------------------------------------------
    

    mpz_t y;
    mpz_init(y);
    mpz_t pooltoseive[MAXVECTOR_SIZE];

    //---------------Constructing the factor base---------
    unsigned int *primes;
    primes = malloc((nprimes)*sizeof(unsigned int));
    findFirstPrimesVerifieQuadraticResidue2(n,primes,nprimes);
    
    unsigned int **primesresidue;
    primesresidue = malloc((nprimes)*sizeof(unsigned int*));
    for (i=1; i<nprimes; i++) {
          mpz_set_ui(tmp,primes[i]);
          findQuadraticResidue(n,tmp,&r1,&r2);
          primesresidue[i]=malloc(2*sizeof(unsigned int));
	  primesresidue[i][0]=mpz_get_ui(r1);
          primesresidue[i][1]=mpz_get_ui(r2);
    }
    
    mpz_t invAprimes[nprimes];
    for (j=1; j<nprimes; j++)
    {
	mpz_init(invAprimes[j]);            
    }
    if (verbose) printf("\nprocess %d : Had generated successfully Factor Base with %ld primes ...\n",rank,nprimes);
    //----------------------------------------------------
    unsigned long columnsize= nprimes/8+1;
    unsigned char **bitarray;
    bitarray = malloc((MAXVECTOR_SIZE)*sizeof(char*));
    for (i=0; i<MAXVECTOR_SIZE; i++) bitarray[i] = malloc((columnsize)*sizeof(char));

    unsigned char **matrixtoreduce;
    matrixtoreduce = malloc((maxalloc)*sizeof(char*));
    for (i=0; i<maxalloc; i++) matrixtoreduce[i] = malloc((columnsize)*sizeof(char));

    /*unsigned char **matrixtoreducecopy;
    matrixtoreducecopy = malloc((maxalloc)*sizeof(char*));
    for (i=0; i<maxalloc; i++) matrixtoreducecopy[i] = malloc((columnsize)*sizeof(char));*/

    //long lineofx[nprimes+maxsmooths];
    long *lineofx;
    lineofx = malloc((maxalloc)*sizeof(long));

    
    //unsigned char lineofxPoly[nprimes+maxsmooths];
    int *lineofxPoly;
    lineofxPoly = malloc((maxalloc)*sizeof(int));
    
    mpz_t polynomeFactors[MAXPOLYNOMES][3];

    unsigned long count=0;
    
    //--------------------LOGARITHM SIEVING-----------------------------
    /*mpz_t x_max;
    double thresh;
    mpz_init(x_max);
    mpz_set(x_max,n);
    mpz_mul_ui(x_max,x_max,2);
    mpz_sqrt(x_max,x_max);
    mpz_mul_ui(x_max,x_max,poolsize);
    mpz_tdiv_q_ui(x_max,x_max,2);
    thresh = (double) logarithm10(x_max);
    thresh=thresh * 0.735;
    double fudge=0;
    long maxprime=(long)(thresh*3);
    for (i=1;i<nprimes;i++) {
        if (primes[i]>=maxprime) break;
        fudge+=log10((double) primes[i]);
        
    }
    fudge=fudge/4;
    thresh-=fudge;*/
    //--------------------LOGARITHM SIEVING-----------------------------
	
//---------------PARALLEL SEGMENT --------------------//

	



//-------------PARALLELE SEGMENT ---------------------//
    timebeginseive = clock();
    
    while (!stops[rank])
    {   
	if (numpoly!=0)		
        {
            /*
            Polynome = (A.X + B)^2 - n
            for sieving while A is a square
            Polynome = A.X^2 + 2.B.X + C
            A = rootnprime * rootnprime
            B = (b + (n - b*b) * mod_inv(b + b, root_A))%A
            B.B-A.C = n <=> C = (B.B-n)/A
            */
            //-----------------Generating polynomes------------------------------------
            
	    //next_prime(rootnprime,&rootnprime);
            //while (!isQuadraticResidue(n,rootnprime)) next_prime(rootnprime,&rootnprime);


	//On a m machine 
	//
	/*le maitre prend le polynome 0
	il va générer les nprime 
	*/



	    
            next_k_prime_QuadraticResidue(n,rootnprime,&rootnprime,nbprocess);
            mpz_mul(A,rootnprime,rootnprime);
            findQuadraticResidue(n,rootnprime,&b,&r2);
            mpz_add(B,b,b);
            mpz_invert(B,B,rootnprime);
            mpz_mul(Btemp,b,b);
            mpz_sub(Btemp,n,Btemp);
            mpz_mul(B,B,Btemp);
            mpz_add(B,B,b);
            mpz_mod(B,B,A);
            mpz_mul(C,B,B);
            mpz_sub(C,C,n);
            mpz_tdiv_q(C,C,A);
            sprintf(s,"(%s . X + %s)^2 - %s, with %s",mpz_get_str (NULL, 10, A),mpz_get_str (NULL, 10, B),mpz_get_str (NULL, 10, n),mpz_get_str (NULL, 10, rootnprime));
            if ((rank==0) && verbose) printf("\n%s\n",s);
            mpz_init(polynomeFactors[numpoly][0]);
            mpz_init(polynomeFactors[numpoly][1]);
            mpz_init(polynomeFactors[numpoly][2]);
            mpz_set(polynomeFactors[numpoly][0],A);
            mpz_set(polynomeFactors[numpoly][1],B);
            mpz_set(polynomeFactors[numpoly][2],C);
            for (j=1; j<nprimes; j++)
            {
		mpz_set_ui(tmp,primes[j]);
                mpz_invert(invAprimes[j],polynomeFactors[numpoly][0],tmp);
                
            }
            //-------------------------------------------------------------------------------
        }
        else
        {
            /*
            First polynome = (X + [racine(n)])^2-n
            */
            sprintf(s,"(X + %s)^2 - %s",mpz_get_str (NULL, 10, rootn),mpz_get_str (NULL, 10, n));

	    //---------------PARALLEL SEGMENT------------------//

		if (rank==0)
		{
			numpoly=0;	
		}
		else
		{		
                        
			numpoly=rank;
			
                        next_k_prime_QuadraticResidue(n,rootnprime,&rootnprime,numpoly);
			
                        mpz_mul(A,rootnprime,rootnprime);
		        findQuadraticResidue(n,rootnprime,&b,&r2);
		        mpz_add(B,b,b);
		        mpz_invert(B,B,rootnprime);
		    	mpz_mul(Btemp,b,b);
		    	mpz_sub(Btemp,n,Btemp);
		    	mpz_mul(B,B,Btemp);
		    	mpz_add(B,B,b);
		    	mpz_mod(B,B,A);
		   	mpz_mul(C,B,B);
		    	mpz_sub(C,C,n);
		    	mpz_tdiv_q(C,C,A);
		    	sprintf(s,"(%s . X + %s)^2 - %s, with %s",mpz_get_str (NULL, 10, A),mpz_get_str (NULL, 10, B),mpz_get_str (NULL, 10, n),mpz_get_str (NULL, 10, rootnprime));
		    	if (rank==0) printf("\n%s\n",s);
		    	mpz_init(polynomeFactors[numpoly][0]);
		    	mpz_init(polynomeFactors[numpoly][1]);
		    	mpz_init(polynomeFactors[numpoly][2]);
		    	mpz_set(polynomeFactors[numpoly][0],A);
		    	mpz_set(polynomeFactors[numpoly][1],B);

		    	mpz_set(polynomeFactors[numpoly][2],C);
		    	for (j=1; j<nprimes; j++)
		   	{
			     mpz_set_ui(tmp,primes[j]);
		             mpz_invert(invAprimes[j],polynomeFactors[numpoly][0],tmp);
		    	}
			
			
		}
	    //---------------PARALLEL SEGMENT------------------//
        }
        unsigned long mulavancer=0;
        long begin=mulavancer*MAXVECTOR_SIZE-poolsize;
        long end=(mulavancer+1)*MAXVECTOR_SIZE-poolsize;
        long i2=0;
        
        
        while (end<poolsize)
        {
            if ((rank==0) && verbose) printf("\nEND = %ld\n",end);
            for (i=0; i<MAXVECTOR_SIZE; i++)
            {
                for (j=0; j<columnsize; j++) bitarray[i][j]=0;
                if (shouldclear==1) mpz_clear(pooltoseive[i]);
                mpz_init(pooltoseive[i]);
                if(numpoly==0)
                {
                    i2 = i+begin;
                    if (i2<0) mpz_sub_ui(y,rootn,-i2);
                    else mpz_add_ui(y,rootn,i2);
                    mpz_pow_ui (y, y, 2);
                    mpz_sub(y,y,n);
                }
                else
                {
                    i2 = i+begin;
                    mpz_set_si(y,i2);
                    mpz_mul(y,y,y);
                    mpz_mul(y,y,polynomeFactors[numpoly][0]);
                    mpz_set_si(tmp,i2);
                    mpz_mul(tmp,tmp,polynomeFactors[numpoly][1]);
                    mpz_mul_ui(tmp,tmp,2);
                    mpz_add(y,y,tmp);
                    mpz_add(y,y,polynomeFactors[numpoly][2]);
                }
                mpz_set(pooltoseive[i],y);
                if (mpz_sgn(pooltoseive[i])==-1)
                {
                    bitarray[i][0]^=(1 << 0);
                    mpz_neg(pooltoseive[i],pooltoseive[i]);
                }
            }

            i=0;
            mpz_t x1,x2;
            mpz_t rem;
            mpz_init(rem);
            unsigned int increment;
            mpz_init(x1);
            mpz_init(x2);

            for (j=1; j<nprimes; j++)
            {
                mpz_set_ui(prime,primes[j]);
                mpz_set_ui(r1,primesresidue[j][0]);
                mpz_set_ui(r2,primesresidue[j][1]);
                if (numpoly>=1) mpz_set(invA,invAprimes[j]);
                if (mpz_cmp_ui(prime,2)==0)
                {
                    increment=primes[j];
                    if (numpoly==0)
                    {
                        mpz_sub(x1,r1,rootn);
                        mpz_mod(x1,x1,prime);
                    }
                    else
                    {
                        /*x1 = int(((mod_root[i] - B) * inv_A)%p)
                        x2 = int(((p - mod_root[i] - B) * inv_A)%p)*/
                        mpz_sub(x1,r1,polynomeFactors[numpoly][1]);
                        mpz_mul(x1,x1,invA);
                        mpz_mod(x1,x1,prime);
                    }
                    i = mpz_get_ui(x1);
                    test2 = increment;
                    test = begin/test2;
                    if (begin!=0) i2=i+(test)*increment;
                    else i2 = i;
                    if (begin<0) i2=i2-increment;
                    if ((i2<begin)) i2+=increment;
                    i=i2-begin;
                    while (i<MAXVECTOR_SIZE)
                    {
                        counting2=mpz_scan1(pooltoseive[i],0);
                        mpz_tdiv_q_2exp(pooltoseive[i],pooltoseive[i],counting2);
                        /*mpz_mod(rem,pooltoseive[i],prime);
                        while (mpz_cmp_ui(rem,0)==0)
                        {
                            mpz_tdiv_q(pooltoseive[i],pooltoseive[i],prime);
                            get_elem2(j,numprime,numbit);
                            bitarray[i][numprime]^=(1 << numbit);
                            mpz_mod(rem,pooltoseive[i],prime);
                        }*/
                        if (mpz_cmpabs_ui(pooltoseive[i],1)==0)
                        {
                            lineofx[count]=i2;
                            lineofxPoly[count]=numpoly;
                            for(j2=0; j2<columnsize; j2++) matrixtoreduce[count][j2]=bitarray[i][j2];
                            count++;
                        }
                        i+=increment;
                        i2+=increment;
                    }
                }
                else
                {
                    if (numpoly==0)
                    {
                        mpz_sub(x1,r1,rootn);
                        mpz_mod(x1,x1,prime);
                        mpz_sub(x2,r2,rootn);
                        mpz_mod(x2,x2,prime);
                    }
                    else
                    {
                        /*x1 = int(((mod_root[i] - B) * inv_A)%p)
                        x2 = int(((p - mod_root[i] - B) * inv_A)%p)*/
                        mpz_sub(x1,r1,polynomeFactors[numpoly][1]);
                        mpz_mul(x1,x1,invA);
                        mpz_mod(x1,x1,prime);
                        mpz_sub(x2,r2,polynomeFactors[numpoly][1]);
                        mpz_mul(x2,x2,invA);
                        mpz_mod(x2,x2,prime);
                    }
		increment=primes[j];
                i = mpz_get_ui(x1);
                i1 = mpz_get_ui(x2);
                test2 = increment;
                test = begin/test2;
                if (begin!=0) i2=i+(test)*increment;
                else i2 = i;
                if (begin<0) i2=i2-increment;
                if ((i2<begin)) i2+=increment;
                i=i2-begin;
                
                if (begin!=0) i2=i1+(test)*increment;
                else i2 = i1;
                if (begin<0) i2=i2-increment;
                if ((i2<begin)) i2+=increment;
                i1=i2-begin;
                 
                if (i<i1) {k=i;inc1=i1-i;}
                else {k=i1;inc1=i-i1;}
                actualincrementer=inc1;
                inc2=increment-inc1; 
                
		while (k<MAXVECTOR_SIZE){
                        i2=k+begin;        
                        mpz_mod(rem,pooltoseive[k],prime);
			//mpz_set_ui(rem,0);
                        while (mpz_cmp_ui(rem,0)==0)
                        {
                           mpz_tdiv_q(pooltoseive[k],pooltoseive[k],prime);
                           get_elem2(j,numprime,numbit);
                           bitarray[k][numprime]^=(1 << numbit);
                           mpz_mod(rem,pooltoseive[k],prime);
                        }
                        if (mpz_cmpabs_ui(pooltoseive[k],1)==0)
                        {
                           lineofx[count]=i2;
                           lineofxPoly[count]=numpoly;
                           for(j2=0; j2<columnsize; j2++) matrixtoreduce[count][j2]=bitarray[k][j2];
                           count++;
                        }      
                    if(actualincrementer==inc1) {k+=inc1;actualincrementer=inc2;}
                    else {k+=inc2;actualincrementer=inc1;}  
                }
                }
            }

            mulavancer=mulavancer+1;
            begin=mulavancer*MAXVECTOR_SIZE-poolsize;
            end=(mulavancer+1)*MAXVECTOR_SIZE-poolsize;
        }



        //----------------Last sieving iteration-----------------------------
        

        unsigned long ending = MAXVECTOR_SIZE - (end-poolsize) ;
        
	for (i=0; i<ending; i++)
        {
            for (j=0; j<columnsize; j++) bitarray[i][j]=0;
            if (shouldclear==1) mpz_clear(pooltoseive[i]);
            shouldclear=1;
            mpz_init(pooltoseive[i]);
            if(numpoly==0)
            {
                i2 = i+begin;
                if (i2<0) mpz_sub_ui(y,rootn,-i2);
                else mpz_add_ui(y,rootn,i2);
                mpz_pow_ui (y, y, 2);
                mpz_sub(y,y,n);
            }
            else
            {
                i2 = i+begin;
                mpz_set_si(y,i2);
                mpz_mul(y,y,y);
                mpz_mul(y,polynomeFactors[numpoly][0],y);
                mpz_set_si(tmp,i2);
                mpz_mul(tmp,tmp,polynomeFactors[numpoly][1]);
                mpz_mul_ui(tmp,tmp,2);
                mpz_add(y,y,tmp);
                mpz_add(y,y,polynomeFactors[numpoly][2]);
            }
            mpz_set(pooltoseive[i],y);
            if (mpz_sgn(pooltoseive[i])==-1)
            {
                bitarray[i][0]^=(1 << 0);
                mpz_neg(pooltoseive[i],pooltoseive[i]);
            }
        }
        i=0;
        mpz_t x1,x2;
        mpz_t rem;
        mpz_init(rem);
        long increment;
        mpz_init(x1);
        mpz_init(x2);
        for (j=1; j<nprimes; j++)
        {
            mpz_set_ui(prime,primes[j]);
            mpz_set_ui(r1,primesresidue[j][0]);
            mpz_set_ui(r2,primesresidue[j][1]);
            if (numpoly>=1) mpz_set(invA,invAprimes[j]);
            if (mpz_cmp_ui(prime,2)==0)
            {
                increment=primes[j];
                if (numpoly==0)
                {
                    mpz_sub(x1,r1,rootn);
                    mpz_mod(x1,x1,prime);
                }
                else
                {
                    mpz_sub(x1,r1,polynomeFactors[numpoly][1]);
                    mpz_mul(x1,x1,invA);
                    mpz_mod(x1,x1,prime);
                }
                i = mpz_get_ui(x1);
                test2 = increment;
                test = begin/test2;
                if (begin!=0) i2=i+(test)*increment;
                else i2 = i;
                if (begin<0) i2=i2-increment;
                if ((i2<begin)) i2+=increment;
                i=i2-begin;

                while (i<ending)
                {
                    
                    /*mpz_set_ui(rem,0);
                    while (mpz_cmp_ui(rem,0)==0)
                    {
                        mpz_tdiv_q(pooltoseive[i],pooltoseive[i],prime);
                        get_elem2(j,numprime,numbit);
                        bitarray[i][numprime]^=(1 << numbit);
                        mpz_mod(rem,pooltoseive[i],prime);
                    }*/
                    counting2=mpz_scan1(pooltoseive[i],0);
                    mpz_tdiv_q_2exp(pooltoseive[i],pooltoseive[i],counting2);
                    if (mpz_cmpabs_ui(pooltoseive[i],1)==0)
                    {
                        lineofx[count]=i2;
                        lineofxPoly[count]=numpoly;
                        for(j2=0; j2<columnsize; j2++) matrixtoreduce[count][j2]=bitarray[i][j2];
                        count++;
                    }
                    i+=increment;
                    i2+=increment;
                }
            }
            else
            {
                if (numpoly==0)
                {
                    mpz_sub(x1,r1,rootn);
                    mpz_mod(x1,x1,prime);
                    mpz_sub(x2,r2,rootn);
                    mpz_mod(x2,x2,prime);
                }
                else
                {
                    mpz_sub(x1,r1,polynomeFactors[numpoly][1]);
                    mpz_mul(x1,x1,invA);
                    mpz_mod(x1,x1,prime);
                    mpz_sub(x2,r2,polynomeFactors[numpoly][1]);
                    mpz_mul(x2,x2,invA);
                    mpz_mod(x2,x2,prime);

                }

                
		increment=primes[j];
                i = mpz_get_ui(x1);
                i1 = mpz_get_ui(x2);
                test2 = increment;
                test = begin/test2;
                if (begin!=0) i2=i+(test)*increment;
                else i2 = i;
                if (begin<0) i2=i2-increment;
                if ((i2<begin)) i2+=increment;
                i=i2-begin;
                
                actualincrementer;
                if (begin!=0) i2=i1+(test)*increment;
                else i2 = i1;
                if (begin<0) i2=i2-increment;
                if ((i2<begin)) i2+=increment;
                i1=i2-begin;
                 
                if (i<i1) {k=i;inc1=i1-i;}
                else {k=i1;inc1=i-i1;}
                actualincrementer=inc1;
                inc2=increment-inc1; 
                
		while (k<ending){
                       
                        i2=k+begin;        
                        mpz_mod(rem,pooltoseive[k],prime);
			while (mpz_cmp_ui(rem,0)==0)
                        {
                           mpz_tdiv_q(pooltoseive[k],pooltoseive[k],prime);
                           get_elem2(j,numprime,numbit);
                           bitarray[k][numprime]^=(1 << numbit);
                           mpz_mod(rem,pooltoseive[k],prime);
                        }
                        if (mpz_cmpabs_ui(pooltoseive[k],1)==0)
                        {
                           lineofx[count]=i2;
                           lineofxPoly[count]=numpoly; 
                           for(j2=0; j2<columnsize; j2++) matrixtoreduce[count][j2]=bitarray[k][j2];
                           count++;
                        }      
                     
                    if(actualincrementer==inc1) {k+=inc1;actualincrementer=inc2;}
                    else {k+=inc2;actualincrementer=inc1;}  
		    
                }
            }
        }
	//----------------Last sieving iteration-----------------------------
        int or=0;
        long wantedx;
        long xpoly;
        mpz_t X;
        mpz_t divider1,divider2;
        mpz_init(divider1);
        mpz_init(divider2);
        mpz_init(X);
        //--------------Finding Perfect Square(Null Vector)------------------
        /*for (i=0; i<count; i++)
        {
            or=0;
            if (!(or=isNullVector(matrixtoreduce[i],columnsize))) continue;
            wantedx=lineofx[i];
            xpoly=lineofxPoly[i];
            printf("\nFound a Null Vector\n");
            if (xpoly==0)
            {
                if (wantedx<0) mpz_sub_ui(y,rootn,-wantedx);
                else mpz_add_ui(y,rootn,wantedx);
                mpz_pow_ui (y, y, 2);
                mpz_sub(y,y,n);
                if (wantedx<0) mpz_sub_ui(X,rootn,-wantedx);
                else mpz_add_ui(X,rootn,wantedx);
                if (mpz_sgn(y)==-1)
                {
                    mpz_neg(y,y);
                    mpz_sqrt(y,y);
                    mpz_neg(y,y);
                }
                else mpz_sqrt(y,y);

            }
            else
            {
                mpz_set_si(y,wantedx);
                mpz_mul(y,y,y);
                mpz_mul(y,polynomeFactors[xpoly][0],y);
                mpz_set_si(tmp,wantedx);
                mpz_mul(tmp,tmp,polynomeFactors[xpoly][1]);
                mpz_mul_ui(tmp,tmp,2);
                mpz_add(y,y,tmp);
                mpz_add(y,y,polynomeFactors[xpoly][2]);
                mpz_mul(y,y,polynomeFactors[xpoly][0]);
                if (mpz_sgn(y)==-1)
                {
                    mpz_neg(y,y);
                    mpz_sqrt(y,y);
                    mpz_neg(y,y);
                }
                else mpz_sqrt(y,y);
                mpz_set_si(X,wantedx);
                mpz_mul(X,X,polynomeFactors[xpoly][0]);
                mpz_add(X,X,polynomeFactors[xpoly][1]);
            }

            mpz_sub(divider1,X,y);
            mpz_gcd(divider2,n,divider1);
            or=0;
            if (!(mpz_cmp_ui(divider2,1)==0 || mpz_cmp(divider2,n)==0))
            {
                mpz_tdiv_q(divider1,n,divider2);
                mpz_set(*p,divider1);
                mpz_set(*q,divider2);
                printf("\nFound a solution with a Null Vector,\nPolynome used = %s\n",s);
                return;
            }
            /*
             if a line exist which refer to a prime square
             else we need to do gaussian reduction to find a combinaison of
             vectors which yield to null vector (are exponent of primes are
             there for even ----> the multiplication of ys which corresponds
             to this vectors is a square
            
        }*/
	//Every secondary process : if he found maxsmooths/k (k = MPI_COMM_SIZE) send lineofx, lineofxpoly, Matrixtoreduce
        //If master : concatene lineofxs lineofxspoly, Matrixtoreduce Generate f polynomes (f = maxordre) 
        //--------------------------Start reducing the matrix if we have more than nprimes Vectors---------
        


	
	if(count > ((maxsmooths/nbprocess)+1)){
		unsigned int *counts;
		int *displ,*maxs;
		if (rank==0)
		{
			counts = malloc(sizeof(unsigned int)*nbprocess);
			displ = malloc(sizeof(int)*nbprocess);
				
			maxs = malloc(sizeof(int)*nbprocess);		
		}
		int max = lineofxPoly[count-1];
		//printf("I'm process %d ill send max: %d\n",rank,max);		
		MPI_Gather(&max, 1, MPI_INT, maxs, 1, MPI_INT, 0,MPI_COMM_WORLD);
 		//printf("I'm process %d i sent max : %d\n",rank,max);
		MPI_Gather((int*) (&count),1,MPI_INT,counts,1,MPI_INT,0,MPI_COMM_WORLD);
		
				
		

		if (rank==0)
		{
			displ[0]=0;
			int uv=1;
			for (uv=1;uv<nbprocess;uv++)
			{
				displ[uv] = displ[uv-1] + counts[uv-1]; //quandfini//
			}
			
		}
		
		MPI_Gatherv(lineofx,count,MPI_LONG,lineofx,counts,displ,MPI_LONG,0,MPI_COMM_WORLD);
		if (rank!=0)
		{
			int tv=0;
			for (tv=0;tv<count;tv++)
			{
				MPI_Send(matrixtoreduce[tv],columnsize,MPI_BYTE,0,0,MPI_COMM_WORLD);
			}
		}
		else
		{
			int matrixtoreduceindex = count;
			int currRank = 0;
			int currIndexInRank = 0;
			for (currRank=1;currRank<nbprocess;currRank++)
			{	
				for (currIndexInRank=0;currIndexInRank<counts[currRank];currIndexInRank++)
				{
					
					MPI_Recv(matrixtoreduce[matrixtoreduceindex],columnsize,MPI_BYTE,currRank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					matrixtoreduceindex++;
					
				}
			}
		}
				
			
		MPI_Gatherv(lineofxPoly,count,MPI_INT,lineofxPoly,counts,displ,MPI_INT,0,MPI_COMM_WORLD);
		if (rank==0)
		{
			count = displ[nbprocess-1] + counts[nbprocess-1];
			max = -1;
			int iv= 0;
			for (iv = 0 ; iv<nbprocess; iv++)
			{
				if (maxs[iv]>max)max=maxs[iv];
						
			}
		}		
		stops[rank] = 1;
		if(rank==0){
			mpz_set(rootnprime,rootnprimebegin);
  			unsigned int generatepoly=1;
			while (generatepoly<=max){
				
				next_k_prime_QuadraticResidue(n,rootnprime,&rootnprime,1);
				mpz_mul(A,rootnprime,rootnprime);
		        	findQuadraticResidue(n,rootnprime,&b,&r2);
		        	mpz_add(B,b,b);
		       	 	mpz_invert(B,B,rootnprime);
		    		mpz_mul(Btemp,b,b);
		    		mpz_sub(Btemp,n,Btemp);
		    		mpz_mul(B,B,Btemp);
		   		mpz_add(B,B,b);
		    		mpz_mod(B,B,A);
		    		mpz_mul(C,B,B);
		    		mpz_sub(C,C,n);
		    		mpz_tdiv_q(C,C,A);
		    		mpz_init(polynomeFactors[generatepoly][0]);
		    		mpz_init(polynomeFactors[generatepoly][1]);
		    		mpz_init(polynomeFactors[generatepoly][2]);
		    		mpz_set(polynomeFactors[generatepoly][0],A);
		    		mpz_set(polynomeFactors[generatepoly][1],B);
		    		mpz_set(polynomeFactors[generatepoly][2],C);
                                //printf("\ngenerated=%d,poly= %s\n",generatepoly,mpz_get_str(NULL,10,polynomeFactors[generatepoly][0]));
				generatepoly++;
			}
			
		}
		/*int MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                void *recvbuf, const int *recvcounts, const int *displs,
                MPI_Datatype recvtype, int root, MPI_Comm comm)*/
	
        }
	if((rank==0) && (!or) && (count > maxsmooths))
        {
            timeendseive = clock();
            
            
            
            
            /*for (i=0;i<count;i++){
            	for (j=0;j<columnsize;j++){
                 matrixtoreducecopy[i][j]=matrixtoreduce[i][j];
               }
            }*/
		
            int trouve=0;
            
            unsigned long f;
            unsigned long marked[nprimes];
            unsigned long markedindex = 0;
            unsigned char testbit=0;
            unsigned char res;
            //unsigned char permuteline[columnsize];
	    unsigned char *permuteline;
            long permute;
            int permutepoly;
	    k=0;
            for (j=0; j<nprimes; j++)
            {
                //if (verbose) printf("\nActual Prime = %ld,%d\n",j,primes[j]);
                if (verbose) printf("\nPercentage of reduction= %.2f\n", 100 * ((float) j/ (float) (nprimes-1)));
                
                k=markedindex;
                get_elem2(j,numprime,numbit);
                while ((((matrixtoreduce[k][numprime] & (1 << numbit)) >> numbit)!=1) && (k<count)) k++;
                if (k<count)
                {
		    /*for (f=0; f<columnsize; f++){
		    	permuteline[f] = matrixtoreducecopy[k][f];
			matrixtoreducecopy[k][f] = matrixtoreducecopy[markedindex][f];
			matrixtoreducecopy[markedindex][f] = permuteline[f];
		    }*/
                    permuteline = matrixtoreduce[k];
                    matrixtoreduce[k] = matrixtoreduce[markedindex];
                    matrixtoreduce[markedindex] = permuteline;
                    permuteline = NULL;
                    
                    permute = lineofx[k];
                    permutepoly = lineofxPoly[k];
                    lineofx[k] = lineofx[markedindex];
                    lineofxPoly[k] = lineofxPoly[markedindex];
                    lineofx[markedindex] = permute; 
                    lineofxPoly[markedindex] = permutepoly;
                    k = markedindex;
                    marked[markedindex]=k;
		    markedindex++;
                    for(i=0; i<nprimes; i++)
                    {
                        if (i!=j)
                        {
                            get_elem2(i,numprime,numbit);
                            
                            if (((matrixtoreduce[k][numprime] & (1 << numbit)) >> numbit)==1)
                            {
                                get_elem2(j,numprime2,numbit2);
                                for (f=(markedindex-1); f<count; f++)
                                {
                                    
                                    matrixtoreduce[f][numprime]^=(((matrixtoreduce[f][numprime2] & (1 << numbit2)) >> numbit2) << numbit);
                                }
                            }
                        }
                    }
                }
            }
            if (verbose) printf("\n Finished with Fast Gauss Reduction\n");
            long basevector[nprimes+1];
            for (i=0; i<markedindex; i++)
            {
                j=0;
                get_elem2(j,numprime,numbit);
                testbit =  (matrixtoreduce[marked[i]][numprime] & (1 << numbit)) >> numbit;
                while (testbit!=1)
                {
                    j++;
                    get_elem2(j,numprime,numbit);
                    testbit =  (matrixtoreduce[marked[i]][numprime] & (1 << numbit)) >> numbit;
                }
                //basevector[i]=j;
		basevector[j]=marked[i];
            }

            long horsbase[count];
            long indexhorsbase=0;
            int boolean=0;
            i=0;
            while (i<count)
            {
                boolean=1;
                for (j=0; j<markedindex; j++)
                {
                    if (i==marked[j])
                    {
                        boolean=0;
                        break;
                    }
                }

                if (boolean==1)
                {
                 
                    horsbase[indexhorsbase]=i;
                    indexhorsbase++;
                }
                i=i+1;
            }
            unsigned long indexselected=0;
            //printf("\nindexhorsbase = %ld\n",indexhorsbase);
            while (!trouve && indexselected<indexhorsbase)
            {
                long select = horsbase[indexselected];
                long selectedx[count];
                unsigned long selectedxindex=0;
                selectedx[selectedxindex]=select;
                selectedxindex=selectedxindex+1;
                for (j=0; j<nprimes; j++)
                {
                    get_elem2(j,numprime,numbit);
                    testbit =  (matrixtoreduce[select][numprime] & (1 << numbit)) >> numbit;
                    if(testbit==1)
                    {
                        f=0;
                        //while(basevector[f]!=j) f++;
                        selectedx[selectedxindex]=basevector[j];
                        selectedxindex=selectedxindex+1;
                    }
                }

                mpz_t modX,modY;
                mpz_init(modX);
                mpz_init(modY);
                mpz_t multiplyofX,multiplyofY;
                mpz_init(multiplyofX);
                mpz_init(multiplyofY);
                mpz_set_ui(multiplyofX,1);
                mpz_set_ui(multiplyofY,1);
                for (i=0; i<selectedxindex; i++)
                {
                    xpoly=lineofxPoly[selectedx[i]];
                    selectedx[i]=lineofx[selectedx[i]];
                    if (xpoly==0)
                    {
                        if (selectedx[i]<0) mpz_sub_ui(X,rootn,-selectedx[i]);
                        else mpz_add_ui(X,rootn,selectedx[i]);
                        mpz_pow_ui (y, X, 2);
                        mpz_sub(y,y,n);
                    }
                    else
                    {
                        mpz_set_si(X,selectedx[i]);
                        mpz_mul(X,X,polynomeFactors[xpoly][0]);
                        mpz_add(X,X,polynomeFactors[xpoly][1]);
                        mpz_set_si(y,selectedx[i]);
                        mpz_mul(y,y,y);
                        mpz_mul(y,y,polynomeFactors[xpoly][0]);
                        mpz_mul(y,y,polynomeFactors[xpoly][0]);
                        mpz_set_si(tmp,selectedx[i]);
                        mpz_mul(tmp,tmp,polynomeFactors[xpoly][1]);
                        mpz_mul_ui(tmp,tmp,2);
                        mpz_mul(tmp,tmp,polynomeFactors[xpoly][0]);
                        mpz_add(y,y,tmp);
                        mpz_mul(tmp,polynomeFactors[xpoly][0],polynomeFactors[xpoly][2]);
                        mpz_add(y,y,tmp);
                    }

                    mpz_mul(modX,X,X);
                    mpz_mod(modX,modX,n);
                    mpz_mod(modY,y,n);
                    mpz_mul(multiplyofX,multiplyofX,X);
                    mpz_mul(multiplyofY,multiplyofY,y);
                    mpz_pow_ui(modX,multiplyofX,2);
                    mpz_mod(modX,modX,n);
                    mpz_mod(modY,multiplyofY,n);
                }
                mpz_mul(modX,multiplyofX,multiplyofX);
                mpz_mod(modX,modX,n);
                mpz_mod(modY,multiplyofY,n);
                if (mpz_sgn(multiplyofY)==-1)
                {
                    mpz_neg(multiplyofY,multiplyofY);
                    mpz_sqrt(y,multiplyofY);
                    mpz_neg(y,y);
                }
                else mpz_sqrt(y,multiplyofY);

                mpz_mul(modX,multiplyofX,multiplyofX);
                mpz_mod(modX,modX,n);
                mpz_mul(modY,y,y);
                mpz_mod(modY,modY,n);
                mpz_sub(divider1,multiplyofX,y);
                mpz_gcd(divider2,n,divider1);
                if (!(mpz_cmp_ui(divider2,1)==0 || mpz_cmp(divider2,n)==0))
                {
                    mpz_tdiv_q(divider1,n,divider2);
                    mpz_set(*p,divider1);
                    mpz_set(*q,divider2);
                    trouve=1;
                    if (verbose) printf("\nPOLYNOME = %s,count = %ld\n",s,count);
                    
		    return;
                    break;
                }
                else
                {
                    indexselected++;
                }
            }
        }
        //printf("\nUsed Polynome : %s\n",s);
        if ((rank==0) && verbose) printf("\nPercentage of sieving : %.2f , NumSmooths = %ld, numpoly = %ld \n", ((double) count/((maxsmooths/nbprocess))) * 100, count,numpoly);
        numpoly = numpoly+nbprocess;
        
        //getchar();
    }
    if (verbose) printf("\nprocess : %d Finished !\n",rank);
}

void factorizeSelfInitializing(mpz_t n,mpz_t *p,mpz_t *q){
    /*long nprimes = (long) logarithm10(n);
    nprimes=5*nprimes*nprimes;
    long poolsize=nprimes*60;
    printf("\nNprimes = %ld,Poolzsize = %ld\n",nprimes,poolsize);*/
    long poolsize;
    long nprimes;
    unsigned int sizebase = mpz_sizeinbase (n,10);
    if (sizebase<=5) {
	nprimes = 20;
        poolsize = 20;  
    }
    else if (sizebase<=10) {
	nprimes = 50;
        poolsize = 500; 
    }
    else if (sizebase<=15) {
	nprimes = 50;
        poolsize = 500; 
    }
    else if (sizebase<=40) {
	nprimes = sizebase * sizebase/ 2;
        poolsize = nprimes * 10;
    }
    /*else if (sizebase<=20) {
	nprimes = 50;
        poolsize = 2500; 
    }
    else if (sizebase<=25) {
	nprimes = 100;
        poolsize = 3000; 
    }
    else if (sizebase<=30) {
	nprimes = 400;
        poolsize = 4000; 
    }
    else if (sizebase<=35) {
	nprimes = 800;
        poolsize = 8000; 
    }*/
    else if (sizebase<=46) {
	nprimes = 167*sizebase-5666;
        poolsize = nprimes*10; 
    }
    else if (sizebase<=46) {
	nprimes = 2000;
        poolsize = 20000; 
    }
    else if (sizebase<=50) {
	nprimes = 4500;
        poolsize = 50000; 
    }
    else if (sizebase<=55) {
	nprimes = 5500;
        poolsize = 100000; 
    }
    else if (sizebase<=60) {
	nprimes = 6500;
        poolsize = 150000; 
    }
    else if (sizebase<=65) {
	nprimes = 7200;
        poolsize = 200000; 
    }
    else if (sizebase<=70) {
	nprimes = 8500;
        poolsize = 250000; 
    }
    else {
    	nprimes = 10000;
        poolsize = 300000; 
    }


    factorize(n,p,q,nprimes,poolsize);
}

int main(int argc, char *argv[])
{

    char argument[100]="";
    unsigned int start=0;
    sprintf(argument,"%s",argv[1]);
    if ((argument[0]=='-') && (argument[1]=='v')) {
	verbose=1;
        start=1;
        argc--;
    }
    else verbose=0;
     
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &nbprocess);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    
    stops = malloc(sizeof(int)*nbprocess);
    int kt=0;
    for (kt=0;kt<nbprocess;kt++) stops[kt]=0;
    mpz_t a,z,b2,c2;
    mpz_init_set_str (a, argv[start+1], 10);
    if (argc>2) {
        
	mpz_init_set_str (b2, argv[start+2], 10);
    	mpz_init_set_str (c2, argv[start+3], 10);
    }
    else {
	mpz_init(b2);
	mpz_init(c2);
        mpz_set_ui(b2,0);
        mpz_set_ui(c2,0);
    }
    
    long nprimes;
    long poolsize;
    nprimes=mpz_get_ui(b2);
    poolsize=mpz_get_ui(c2);
    mpz_t p,q;
    mpz_init(p);
    mpz_init(q);
    mpz_t b,res;
    mpz_init(res);
    if(rank==0) {
        /*printf("\nNumber Process = %d",nbprocess);
     	printf("\nNumber to factorize : %s",mpz_get_str (NULL, 10,a));
    	printf("\nCompose de : %ld bits et de %ld digits",mpz_sizeinbase (a, 2),mpz_sizeinbase (a, 10));*/
    }
    double begin,end;

    
    begin = clock();
    if (argc>2) factorize(a,&p,&q,nprimes,poolsize);
    else factorizeSelfInitializing(a,&p,&q);
    if (rank==0) printf("%s",mpz_get_str (NULL, 10,p));
    if (rank==0) printf("\n%s",mpz_get_str (NULL, 10,q));
    end = clock();
    if (rank==0) printf("\n%f\n",(end-begin)/CLOCKS_PER_SEC);
    
    
    MPI_Finalize();
}
