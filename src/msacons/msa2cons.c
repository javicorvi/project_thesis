/* Morten Nielsen mniel@cbs.dtu.dk. March 2009 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"

int     	p_cluster;
float		p_seqid;
WORD		p_alphabet;
int		p_verbose;
int             p_fasta;
int		p_cons;
float   p_wlowcount;
FILENAME        p_blrealmat;
FILENAME        p_reducedalphabet;
int		p_red;
int		p_printffreq;

PARAM   param[] = {
        "-c", VSWITCH   p_cluster, "Perform sequence clustering", "0",
	"-i", VFLOAT	p_seqid, "Percent id for clustering", "0.62",
	"-a",	VWORD	p_alphabet, "Aminoacid alphabet", "ARNDCQEGHILKMFPSTWYV",
	"-v", VSWITCH	p_verbose, "Verbose mode", "0",
	"-fsa", VSWITCH p_fasta, "Input is fasta format", "0",          
	"-cons", VINT	p_cons, "Type of conservation measure [0] Shannon, [1] Kullback [2] Max freq", "0",
	"-wlc", VFLOAT	p_wlowcount, "Weight on pseudo count (prior)", "200",
	"-blf", VFNAME	p_blrealmat, "Real value blosum matrix filename", "/home/mniel/data/matrices/blosum62.freq_rownorm",
	"-ared", VFNAME	p_reducedalphabet, "File with reduced alphabet", "/home/mniel/data/alphabet.logo", 
	"-r", VSWITCH	p_red, "Used reduced amino acids alphabet", "0",
	"-pf", VSWITCH 	p_printffreq, "Print amino acid frequencies", "0",
        0
};

float   *set_bgfreq_aa()

{
        float   *aabgfrec;

        aabgfrec = fvector(0, 19 );

        aabgfrec[0]  = 0.074;
        aabgfrec[1]  = 0.052;
        aabgfrec[2]  = 0.045;
        aabgfrec[3]  = 0.054;
        aabgfrec[4]  = 0.025;
        aabgfrec[5]  = 0.034;
        aabgfrec[6]  = 0.054;
        aabgfrec[7]  = 0.074;
        aabgfrec[8]  = 0.026;
        aabgfrec[9]  = 0.068;
        aabgfrec[10] = 0.099;
        aabgfrec[11] = 0.058;
        aabgfrec[12] = 0.025;
        aabgfrec[13] = 0.047;
        aabgfrec[14] = 0.039;
        aabgfrec[15] = 0.057;
        aabgfrec[16] = 0.051;
        aabgfrec[17] = 0.013;
        aabgfrec[18] = 0.032;
        aabgfrec[19] = 0.073;

        return( aabgfrec );
}

typedef struct klist {
	struct klist *next;
	char	*seq;
	int	n;
	PIRLIST	*sl;
	struct klist *l;
} KLIST;

KLIST	*klist_alloc2()

{
	KLIST	*n;

	if ( ( n = ( KLIST * ) malloc ( sizeof( KLIST ))) != NULL ) {
		n->next = NULL;
		n->seq = NULL;
		n->l = NULL;
		n->sl = NULL;
		n->n = 0;
	}

	return( n );
}

KLIST	*homolog2( char *pep, KLIST *klist, float thr, int l )

{
	KLIST	*kl;
	int	nid;
	int	i, neq;

	neq = (int)(l * thr);

	for ( kl=klist; kl; kl=kl->next ) {

		for ( i=0, nid=0; i<l && nid <= neq ; i++ )
			nid += ( pep[i] != '-' && pep[i] == (kl->seq)[i] );

		if ( nid > neq )
			return( kl );
	}

	return( NULL );
}

float	cal_w_cluster_mi( PIRLIST *slist, float thr )

{
	PIRLIST	*sl, *sl1;
	int	i, o, len;
	float	w, nseq;
	KLIST	*klist, *kl, *h, *lnext, *klnext, *l;
	char	*pep;
	
	len = slist->len;

	klist = NULL;

	printf( "# Cluster ID %f\n", thr );

	nseq = 0.0;

	for ( i=0,sl=slist; sl; sl=sl->next, i++ ) {

		pep = sl->seq;

		kl = klist_alloc2();

		if ( ! kl ) {
			printf( "Error. Cannot allocate klist element\n" );
			exit( 1 );
		}

		kl->seq = pep;
		kl->sl = sl;

		if ( ( h = homolog2( pep, klist, thr, len ) ) != NULL ) {

			kl->l = h->l;
			h->l = kl;
			h->n++;
		}
		else {

			nseq += 1.0;

			kl->n = 1;

			if ( klist ) {
				kl->next = klist;
				klist = kl;
			}
			else {
				klist = kl;
			}
		}

	}

	for ( kl = klist; kl; kl=kl->next ) {
		w = 1.0/kl->n;
		sl1 = kl->sl;
		sl1->w = w;

		for ( l=kl->l; l; l=l->l ) {
			sl1 = l->sl;
			sl1->w = w;

		}
	}

	for ( kl = klist; kl; kl=klnext ) {
		klnext = kl->next;
		for ( l=kl->l; l; l=lnext ) {
			lnext = l->l;

			free(( char * ) l );

		}

		free(( char * ) kl );

	}

	return( nseq );
}

int main(int argc, char * argv[])

{

	PIRLIST	*pirlist, *pir;
	int	ix;
	int	len;
	int 	i, j;
	char 	*alphabet;
	float	nseq;
	float	*nseq_vec;
	int	alen, n;
	int	maxlen;
	FSALIST	*fsalist;
	float	**aa, *nn, cons;
	int	nprot;
	float	*bg;
	float	**blf;
	char	*bl_alphabet;
	int	alen_red;
	char	**alphabet_red;
	float	**aa_red, *bg_red;
	int	found;
	LINELIST	*linelist, *ln;
	int	*translate;
	int	nc;

	pparse( &argc, &argv, param, 1, "pir_file" );

        if ( p_fasta ) {
                fsalist = fsalist_read( argv[1] );

                pirlist = fsalist2pirlist( fsalist );

                fsalist_free( fsalist );
        }
        else
		pirlist = pirlist_read( argv[1] );

	if ( ! pirlist ) {
		printf( "Error. No PIR elements read from file %s\n", argv[1] );
		exit( 1 );
	}

	blf = read_realblosum( p_blrealmat, &bl_alphabet );

	if ( strncmp( bl_alphabet, PROFILE_ORDER, MIN( strlen( bl_alphabet ), strlen( PROFILE_ORDER ) ) - 1 ) != 0 ) {
                printf( "Error. Bl alphabet %s not equal to %s\n",
                        bl_alphabet, PROFILE_ORDER );
                exit( 1 );
        }

	alen = strlen( p_alphabet );
	alphabet = cvector( 0, alen );
	strcpy( alphabet, p_alphabet );

	printf( "# Amino acid alphabet %s\n", alphabet );

	for ( maxlen=0,nprot=0,pir = pirlist; pir; pir=pir->next, nprot++ ) {
		if ( pir->len > maxlen )
			maxlen = pir->len;
	}

	printf( "# Nprot %i. Maxlen %i\n", nprot, maxlen );

   	len = pirlist->len;

	printf( "# Len %i\n", len );

	if ( p_cluster ) 
        	nseq = cal_w_cluster_mi( pirlist, p_seqid );
	else
		nseq = nprot;

	printf( "# Nseq %f\n", nseq );

	aa = fmatrix( 0, len-1, 0, alen-1 );
	nn = fvector( 0, len-1 );

	nseq_vec = fvector( 0, len-1 );

	set_NAA( alen );

	for ( i=0; i<pirlist->len; i++ ) {

		for ( pir=pirlist; pir; pir=pir->next ) {

			ix = strpos( alphabet, pir->seq[i] );

			if ( ix >=0 )
				nseq_vec[i] += pir->w;

		}

	}

	for ( pir=pirlist; pir; pir=pir->next ) {

		for ( i=0; i<pir->len; i++ ) {

			ix = strpos( alphabet, pir->seq[i] );

			if ( ix >=0 ) {
				aa[i][ix] += pir->w;

				nn[i] += pir->w;
			}
		}
	}

	if ( p_cons == 1 )
		bg = set_bgfreq_aa();
	else
		bg = NULL;

	for ( i=0; i<len; i++ ) {

		for ( j=0; j<alen; j++ )
			aa[i][j] /= ( nn[i] > 0 ? nn[i] : 1.0 );

	}

	if ( p_wlowcount > 0.0 ) {
		lowcount_correct_nseq_vec( aa, nseq_vec, len, blf, p_wlowcount );
	}

	if ( p_red ) {

		linelist = linelist_read( p_reducedalphabet );

		if ( ! linelist ) {
			printf( "Error. Cannot read file with reduced alphabet from %s. Exit\n", p_reducedalphabet );
			exit( 1 );
		}

		for ( alen_red=0, ln = linelist; ln; ln=ln->next, alen_red++ );

		alphabet_red = cmatrix( 0, alen_red-1, 0, alen );

		for ( i=0, ln = linelist; ln; ln=ln->next, i++ ) {
			if ( ( nc = sscanf( ln->line, "%*s %s", alphabet_red[i] ) ) != 1 ) {
				printf( "Error. Wrong line format for reduced alphabet %s. Exit\n", ln->line );
				exit( 1 );
			}
		}

		printf( "# Reduced alphabet read from %s. Number of classes %i\n", p_reducedalphabet, alen_red );

		for ( i=0; i<alen_red; i++ )
			printf( "# Class %i AA %s\n", i, alphabet_red[i] );

		translate = ivector( 0, alen-1 );

		for ( i=0;i<alen; i++ ) {
			found = 0;
			for ( j=0; ! found && j<alen_red; j++ )
				if ( strpos( alphabet_red[j], alphabet[i] ) >= 0 ) {
					translate[i] = j;
					found = 1;
				}

			if ( ! found ) {
				printf( "Error %c not found in reduced alphabet. Exit\n", alphabet[i] );
				exit( 1 );
			}
		}

		aa_red = fmatrix( 0, len-1, 0, alen_red-1 );

		for ( i=0; i<len; i++ ) 
			for ( j=0; j<alen; j++ ) {
				ix = translate[j];

				aa_red[i][ix] += aa[i][j];
			}

		if ( bg ) {
			bg_red = fvector( 0, alen_red-1 );

			for ( j=0; j<alen; j++ ) {

				ix = translate[j];
				bg_red[ix] += bg[j];
			}
		}
	}
	else {
		aa_red = aa;
		bg_red = bg;
		alen_red = alen;
	}

	if ( p_printffreq ) {
		printf( "# Alphabet           " );

		if ( p_red ) {
			for ( j=0; j<alen_red; j++ )
				printf( " %6s", alphabet_red[j] );
		}
		else
			for ( j=0; j<alen_red; j++ )
				printf( " %6c", alphabet[j] );

		printf( "\n" );
	}	

	for ( i=0; i<len; i++ ) {

		cons = 0.0;

		if ( p_cons == 0 ) { 
			for ( j=0; j<alen_red; j++ ) {

				cons += ( aa_red[i][j] > 0 ? aa_red[i][j]*log(aa_red[i][j]) : 0.0 );
			}

			if ( cons > 0 ) {
				cons += log(alen_red); 

				cons = cons/log(2);
			}
		} 
		else if ( p_cons == 1 ) {

			for ( j=0; j<alen_red; j++ )
				cons += ( aa_red[i][j] > 0 ? aa_red[i][j]*log( aa_red[i][j]/bg_red[j] ) : 0.0 );

			if ( cons < 0 )
				cons = 0.0;

			cons = cons/log(2);

		}
		else {

			cons = fvector_max( aa_red[i], alen_red );

		}

		printf( "%5i %c %6.4f", i, pirlist->seq[i], cons );

		if ( p_printffreq ) {
			printf( " AAfreq" );
			for ( j=0; j<alen_red; j++ )
				printf( " %6.4f", aa_red[i][j] );
		}

		printf( "\n" );
	}

	exit( 0 );
    
}
