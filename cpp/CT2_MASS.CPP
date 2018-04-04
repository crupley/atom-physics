
/*From mcrogan@sgi.com Thu Sep 17 14:57:23 1998
Date: Thu, 17 Sep 1998 13:07:04 -0500 (CDT)
From: Michael Crogan <mcrogan@sgi.com>
To: zmw@physics.purdue.edu
Subject: center_of_mass.c with file i/o support, dynamic memory allocation, and override for matrix size via the command line

/*
 * center_of_mass.c
 *
 * Override matrix size (make sure the input file is large enough first) by
 * running the program with max_x as the first argument and max_y as the
 * second, as in:
 *
 * center_of_mass 180 180
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FALSE 0
#define TRUE 1

#define INPUT_FILENAME "c:\\tc\\image134\\img2237.u01"
//#define INPUT_FILENAME "input"
#define OUTPUT_FILENAME "output"

#define MAX_NUMBER_OF_GROUPS 180

#define MAX_X 122
#define MAX_Y 122

int search_matrix (void);
int traverse_group (int x, int y);

int mass = 0;
int group = 0;

int max_x = MAX_X;
int max_y = MAX_Y;

int **probed;
int **data;

int group_x [MAX_NUMBER_OF_GROUPS];
int group_y [MAX_NUMBER_OF_GROUPS];

int weighted_mass [2] = {0};

int main (int argc, char **argv)
  {
    int i, j;
    FILE *input_fp;
    FILE *output_fp;

    if (argc >= 2)
      {
	max_x = atoi (argv [1]);
	if (max_x <= 0)
	  max_x = 1;
      }

    if (argc >= 3)
      {
	max_y = atoi (argv [2]);
	if (max_y <= 0)
	  max_y = 1;
      }

    if ((input_fp = fopen (INPUT_FILENAME, "r")) == NULL)
      {
	fprintf (stderr,
	  "unable to open input file \"" INPUT_FILENAME "\"\n");
	return 1;
      }

    if ((output_fp = fopen (OUTPUT_FILENAME, "w")) == NULL)
      {
	fprintf (stderr,
	  "unable to open output file \"" OUTPUT_FILENAME "\"\n");
	return 1;
      }

    if ((data = (int **) calloc (max_x, sizeof (int *))) == NULL)
      {
	fprintf (stderr,
	  "calloc (%d, %d) fails for data\n", max_x, sizeof (int *));
	return 1;
      }

    for (i = max_x; i--;)
      if ((data [i] = (int *) calloc (max_y, sizeof (int))) == NULL)
	{
	  fprintf (stderr,
	    "calloc (%d, %d) fails for data at [%d]\n", max_y, sizeof (int), i);
	  return 1;
	}

    if ((probed = (int **) calloc (max_x, sizeof (int *))) == NULL)
      {
	fprintf (stderr,
	  "calloc (%d, %d) fails for probed\n", max_x, sizeof (int *));
	return 1;
      }

    for (i = max_x; i--;)
      if ((probed [i] = (int *) calloc (max_y, sizeof (int))) == NULL)
	{
	  fprintf (stderr,
	    "calloc (%d, %d) fails for probed at [%d]\n", max_y, sizeof (int), i);
	  return 1;
	}

    for (i = max_x; i--;)
      for (j = max_y; j--;)
	if (fscanf (input_fp, "%d", &data [i] [j]) != 1)
	  {
	    fprintf (stderr,
	      "fscanf fails at [%d] [%d]\n", i, j);
	    return 1;
	  }

    fclose (input_fp);

    search_matrix ();

    for (i = max_x; i--;)
		//menset(data [i],(unsigned char)0, max_y *sizeof(int));
		  for(j=max_y; j--;)
		     data[i][j]=0;
    for (i = 0; i < group; i++)
      {
	data [group_x [i]] [group_y [i]] = 1;
      }

    for (i = max_x; i--;)
      {
	for (j = max_y; j--;)
	  fprintf (output_fp, "%d ", data [i] [j]);
	fprintf (output_fp, "\n");
      }

    fclose (output_fp);

    for (i = max_x; i--;)
      free (data [i]);

    free (data);

    return 0;
  }

int search_matrix (void)
  {
    int x;
    int y;

    for (x = max_x; x--;)
      for (y = max_y; y--;)
	if (data [x] [y] && !probed [x] [y])
	  {
	    traverse_group (x, y);
				group_x [group] = weighted_mass [0] / mass;
				//group_x [group] =((double)weighted_mass [0] / mass)+0.5;
				group_y [group++] = weighted_mass [1] / mass;
				//group_y [group++] =((double)weighted_mass [1] / mass)+0.5;

	    mass = 0;
	    weighted_mass [0] = 0;
	    weighted_mass [1] = 0;
	  }

    return 0;
  }

int traverse_group (int x, int y)
  {
    int cache;

    if (! (cache = data [x] [y]) || probed [x] [y])
      {
        return FALSE;
      }

    mass += cache;
    weighted_mass [0] += x * cache;
    weighted_mass [1] += y * cache;
    probed [x] [y] ++;

    if (x)
      traverse_group (x - 1, y);
    if (y)
      traverse_group (x, y - 1);
    if (x != max_x - 1)
      traverse_group (x + 1, y);
    if (y != max_y - 1)
      traverse_group (x, y + 1);

    return TRUE;
  }
