
/*From mcrogan@sgi.com Sun Sep 20 16:47:55 1998
Date: Sat, 19 Sep 1998 19:58:28 -0500 (CDT)
From: Michael Crogan <mcrogan@sgi.com>
To: Zhengmin Wang <zmw@physics.purdue.edu>
Subject: Re: run 118x118

Zhengmin Wang writes:
 >  Hi! Mike:
 > I run your program on 486 in our lab. It can open my real image(180x180)
 > and can create output image with maximum size 118x118. When the size is
 > over 118x118, it shows:
 >
 > 119x119      Calloc(119,2) fails for probed at [1]
 >
 > 120x120            (120,2)                     [5]
 >
 > 125x125            (125,2)                     [25]
 >
 > 128x128            (128,2)                     [36]
 >
 > Can you increase it to 180x180?

There are no inherent limitations in the program.  The limitation is
in the compiler or in the operating system or in the physical RAM
available.  The suggestion is to use "unsigned char" as a data type,
meaning that input data should range from 0 to 255 in magnitude.

The MAX_NUMBER_OF_GROUPS setting could possibly be decreased, but make
sure there are never more groups than MAX_NUMBER_OF_GROUPS.

 > What are 2 in ( ) and what are 1  5  25  36  in  [ ] [ ]......?

The 2 is the size (number of bytes) of an integer, while the
119,120,125,128 are the number of bytes which are being allocated
simultaneously.  This form is chosen because this is the exact
calloc() system call which is being made.  The bracketed [] value
shows which one-dimensional vector caused the memory allocation
failure.  The code happens to number its vectors from max to min (and
not vice versa), so the higher-valued matrices (128x128) only got as
"low" as 36, meaning (128-36)*128 cells in the two-dimensional matrix
were successfully created for the 128x128 case.

 > Thank you for your help!

Sure!  I am having a little trouble with my mechanics (physics 410)
problems.  Can I go over them sometime?

Thanks,

 > Wang

Mike */

/*
 * center_of_mass.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloc.h>

#define FALSE 0
#define TRUE 1

#define INPUT_FILENAME "c:\\tc\\image134\\img2237.u01"
#define OUTPUT_FILENAME "output"

#define MAX_NUMBER_OF_GROUPS 1024

#define MAX_X 167
#define MAX_Y 167

int search_matrix (void);
int traverse_group (int x, int y);

int mass = 0;
int group = 0;

int max_x = MAX_X;
int max_y = MAX_Y;

unsigned char **probed;
unsigned char **data;

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

    if ((data = (unsigned char **) farcalloc (max_x, sizeof (unsigned char *))) == NULL)
      {
        fprintf (stderr,
	  "farcalloc (%d, %d) fails for data\n", max_x, sizeof (unsigned char *));
        return 1;
      }

    for (i = max_x; i--;)
      if ((data [i] = (unsigned char *) farcalloc (max_y, sizeof (unsigned char))) == NULL)
        {
          fprintf (stderr,
	    "farcalloc (%d, %d) fails for data at [%d]\n", max_y, sizeof (unsigned char), i);
          return 1;
	}

    if ((probed = (unsigned char **) farcalloc (max_x, sizeof (unsigned char *))) == NULL)
      {
        fprintf (stderr,
	  "farcalloc (%d, %d) fails for probed\n", max_x, sizeof (unsigned char *));
        return 1;
      }

    for (i = max_x; i--;)
      if ((probed [i] = (unsigned char *) farcalloc (max_y, sizeof (unsigned char))) == NULL)
        {
          fprintf (stderr,
	    "farcalloc (%d, %d) fails for probed at [%d]\n", max_y, sizeof (unsigned char), i);
          return 1;
        }

    {
      int value;
      for (i = max_x; i--;)
        for (j = max_y; j--;)
          {
            if (fscanf (input_fp, "%d", &value) != 1)
              {
                fprintf (stderr,
                  "fscanf fails at [%d] [%d]\n", i, j);
                return 1;
              }
            data [i] [j] = (unsigned char) value;
          }
    }

    fclose (input_fp);

    search_matrix ();

    for (i = max_x; i--;)
      //bzero (data [i], max_y * sizeof (unsigned char));
	for(j=max_y;j--;)
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
            group_x [group] = (weighted_mass [0] + (mass >> 1)) / mass; 
            group_y [group++] = (weighted_mass [1] + (mass >> 1)) / mass; 

            mass = 0;
            weighted_mass [0] = 0;
            weighted_mass [1] = 0;
          }

    return 0;
  }

int traverse_group (int x, int y)
  {
    unsigned char cache;

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
