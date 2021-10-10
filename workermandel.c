/*
 * Polytech Lille 
 * Calcul de l'ensemble de Mandelbrot, Version séquentielle
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>	/* chronometrage */
#include <string.h>     /* pour memset */
#include <math.h>
#include <sys/time.h>
#include <mpi.h>
#include <unistd.h>
#include "rasterfile.h"


/*** 
 * Par Francois-Xavier MOREL (M2 SAR, oct2009): 
 */
/**
 * Convertion entier (4 octets) LINUX en un entier SUN
 * @param i entier à convertir
 * @return entier converti
 */

double my_gettimeofday() {
  struct timeval tmp_time;
  gettimeofday( & tmp_time, NULL);
  return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}

int swap(int i) {
  int init = i;
  int conv;
  unsigned char * o, * d;

  o = ((unsigned char * ) & init) + 3;
  d = (unsigned char * ) & conv;

  * d++ = * o--;
  * d++ = * o--;
  * d++ = * o--;
  * d++ = * o--;

  return conv;
}
unsigned char power_composante(int i, int p) {
  unsigned char o;
  double iD = (double) i;

  iD /= 255.0;
  iD = pow(iD, p);
  iD *= 255;
  o = (unsigned char) iD;
  return o;
}

unsigned char cos_composante(int i, double freq) {
  unsigned char o;
  double iD = (double) i;
  iD = cos(iD / 255.0 * 2 * M_PI * freq);
  iD += 1;
  iD *= 128;

  o = (unsigned char) iD;
  return o;
}

/*** 
 * Choix du coloriage : definir une (et une seule) des constantes
 * ci-dessous :  
 */
//#define ORIGINAL_COLOR
#define COS_COLOR

#ifdef ORIGINAL_COLOR
#define COMPOSANTE_ROUGE(i)((i) / 2)
#define COMPOSANTE_VERT(i)((i) % 190)
#define COMPOSANTE_BLEU(i)(((i) % 120) * 2)
#endif /* #ifdef ORIGINAL_COLOR */
#ifdef COS_COLOR
#define COMPOSANTE_ROUGE(i) cos_composante(i, 13.0)
#define COMPOSANTE_VERT(i) cos_composante(i, 5.0)
#define COMPOSANTE_BLEU(i) cos_composante(i + 10, 7.0)
#endif /* #ifdef COS_COLOR */

/**
 *  Sauvegarde le tableau de données au format rasterfile
 *  8 bits avec une palette de 256 niveaux de gris du blanc (valeur 0)
 *  vers le noir (255)
 *    @param nom Nom de l'image
 *    @param largeur largeur de l'image
 *    @param hauteur hauteur de l'image
 *    @param p pointeur vers tampon contenant l'image
 */

void sauver_rasterfile(char * nom, int largeur, int hauteur, unsigned char * p) {
  FILE * fd;
  struct rasterfile file;
  int i;
  unsigned char o;

  if ((fd = fopen(nom, "w")) == NULL) {
    printf("erreur dans la creation du fichier %s \n", nom);
    exit(1);
  }

  file.ras_magic = swap(RAS_MAGIC);
  file.ras_width = swap(largeur); /* largeur en pixels de l'image */
  file.ras_height = swap(hauteur); /* hauteur en pixels de l'image */
  file.ras_depth = swap(8); /* profondeur de chaque pixel (1, 8 ou 24 )   */
  file.ras_length = swap(largeur * hauteur); /* taille de l'image en nb de bytes		*/
  file.ras_type = swap(RT_STANDARD); /* type de fichier */
  file.ras_maptype = swap(RMT_EQUAL_RGB);
  file.ras_maplength = swap(256 * 3);

  fwrite( & file, sizeof(struct rasterfile), 1, fd);

  /* Palette de couleurs : composante rouge */
  i = 256;
  while (i--) {
    o = COMPOSANTE_ROUGE(i);
    fwrite( & o, sizeof(unsigned char), 1, fd);
  }

  /* Palette de couleurs : composante verte */
  i = 256;
  while (i--) {
    o = COMPOSANTE_VERT(i);
    fwrite( & o, sizeof(unsigned char), 1, fd);
  }

  /* Palette de couleurs : composante bleu */
  i = 256;
  while (i--) {
    o = COMPOSANTE_BLEU(i);
    fwrite( & o, sizeof(unsigned char), 1, fd);
  }

  // pour verifier l'ordre des lignes dans l'image : 
  //fwrite( p, largeur*hauteur/3, sizeof(unsigned char), fd);

  // pour voir la couleur du '0' :
  // memset (p, 0, largeur*hauteur);

  fwrite(p, largeur * hauteur, sizeof(unsigned char), fd);
  fclose(fd);
}



unsigned char xy2color(double a, double b, int prof) {
  double x, y, temp, x2, y2;
  int i;

  x = y = 0.;
  for (i = 0; i < prof; i++) {
    /* garder la valeur précédente de x qui va etre ecrase */
    temp = x;
    /* nouvelles valeurs de x et y */
    x2 = x * x;
    y2 = y * y;
    x = x2 - y2 + a;
    y = 2 * temp * y + b;
    if (x2 + y2 > 4.0) break;
  }
  return (i == prof) ? 255 : (int)((i % 255));
}


#define SIZE_H_N 50
#define NB_LIGNES 50


int main(int argc, char * argv[]) {

    /* Domaine de calcul dans le plan complexe */
    double xmax, ymax;
    double xmin, ymin;
    /* Dimension de l'image */
    int w, h;
    /* Pas d'incrementation */
    double xinc, yinc;
    /* Profondeur d'iteration */
    int prof;
    /* Image resultat */
    unsigned char * ima;
    /* Variables intermediaires */
    int i, j;
    double x, y;
    /* Chronometrage */
    double debut, fin;
    int portion_traitee=0;
    /* debut du chronometrage */
    debut = my_gettimeofday();

    int my_rank; /*rang du processus	*/
    int p; /*nombre de processus*/	
    int tag = 0; /*etiquette du message	*/
    MPI_Status status;

    /* Valeurs par defaut de la fractale */
    xmin = -2;
    ymin = -2;
    xmax = 2;
    ymax = 2;
    w = h = 800;
    prof = 15000;
    /* Recuperation des parametres */
    if (argc > 1) w = atoi(argv[1]);
    if (argc > 2) h = atoi(argv[2]);
    if (argc > 3) xmin = atof(argv[3]);
    if (argc > 4) ymin = atof(argv[4]);
    if (argc > 5) xmax = atof(argv[5]);
    if (argc > 6) ymax = atof(argv[6]);
    if (argc > 7) prof = atoi(argv[7]);
    /* Calcul des pas d'incrementation */
    xinc = (xmax - xmin) / (w - 1);
    yinc = (ymax - ymin) / (h - 1);

    /* affichage parametres pour verificatrion */
    fprintf(stderr, "Domaine: {[%lg,%lg]x[%lg,%lg]}\n", xmin, ymin, xmax, ymax);
    fprintf(stderr, "Increment : %lg %lg\n", xinc, yinc);
    fprintf(stderr, "Prof: %d\n", prof);
    fprintf(stderr, "Dim image: %dx%d\n", w, h);

    MPI_Init( &argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, & p);
    

    if (my_rank == 0){
        /* Allocation memoire du tableau resultat */
        unsigned char *ima = (unsigned char * ) malloc(w*h*sizeof(unsigned char));

        for (int i = 1; (i < p)&&(i<h/NB_LIGNES);i++){
            MPI_Send(MPI_BOTTOM, 0,MPI_UNSIGNED_CHAR,i,portion_traitee,MPI_COMM_WORLD);
            portion_traitee++;
        }
        for (int i=0;i<h/NB_LIGNES;i++){
            unsigned char *cp = (unsigned char *) malloc(w*h*sizeof(unsigned char));
            MPI_Recv(cp,NB_LIGNES*w*sizeof(unsigned char),MPI_UNSIGNED_CHAR,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
            memcpy(&ima[status.MPI_TAG*NB_LIGNES*w],cp,NB_LIGNES*w);
            if (portion_traitee < h/NB_LIGNES){
                MPI_Send(MPI_BOTTOM,0,MPI_UNSIGNED_CHAR,status.MPI_SOURCE,portion_traitee, MPI_COMM_WORLD);
                portion_traitee++;
            } else {
                MPI_Send(MPI_BOTTOM,0,MPI_UNSIGNED_CHAR, status.MPI_SOURCE, h/NB_LIGNES, MPI_COMM_WORLD);
                
            }   
        }
        /* Sauvegarde de la grille dans le fichier resultat "mandel.ras" */
        sauver_rasterfile("mandel.ras", w, h, ima);
        /* fin du chronometrage */
	      fin = my_gettimeofday();
	      fprintf(stderr, "Temps total de calcul : %g sec\n", fin - debut);
    } else {
        unsigned char *tab = (unsigned char * ) malloc(w*NB_LIGNES*sizeof(unsigned char));
        MPI_Recv(&tab, NB_LIGNES*w*sizeof(unsigned char),MPI_UNSIGNED_CHAR,0,MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        int portion = status.MPI_TAG;
        while (portion < h/NB_LIGNES){
            y = ymin + portion*NB_LIGNES*yinc;
            for (int j=0; j<NB_LIGNES*w;j++){
                if (j%w == 0)
                    x = xmin;
                tab[j] = xy2color(x,y,prof);
                x = x + xinc;
                if (j%w == 0)
                    y = y + yinc;
            }
            MPI_Send(tab, NB_LIGNES*w*sizeof(unsigned char),MPI_UNSIGNED_CHAR,0,portion,MPI_COMM_WORLD);
            MPI_Recv(&tab, NB_LIGNES*w*sizeof(unsigned char),MPI_UNSIGNED_CHAR,0,MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            portion = status.MPI_TAG;
            
            
        }
    }
    MPI_Finalize();
    
}