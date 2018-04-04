/*
 From mcrogan@csociety.purdue.edu Mon Sep 28 17:48:15 1998
Date: Sun, 27 Sep 1998 12:37:08 -0500
From: Michael Lee Crogan <mcrogan@csociety.purdue.edu>
To: mcrogan@sgi.com, zmw@physics.purdue.edu
Subject: Integrated code

Wang:

I have tried to integrate Tim's code with mine a bit.  I still think it's
a big priority to check on the change of compilers, but here's what the
integrated code looks like.  The other compiler would use a different
memory allocation method for one.  The biggest reason for switching is
memory.  I sent you a copy of the quick email I composed to Tim.  But
here is the updated code.  Note I have no way to compile this because I
do not have Borland C++.  */

/* Michael Rifani, Aug. 2, 1996                                       */

/* the purpose of this C++ program is to try two externally triggerred */
/* each resulting image in UNDERLAY memory is read by u_get_pixel and  */
/* to test u_get_freg to figure out when one frame is done being read  */
/* and thus the next one can start being read.                         */

/* technical assistant on the Falcon stuff was provided by Ron Green   */
/* of Univision Tech. (tel.: (508) 667-8900, fax: (508) 670-1960),     */
/* and faxes have been more effective in getting their attention       */

/* according to Ron Green:                                             */
/* to compile and link properly, create a Turbo C++ project in which   */
/* the library c:\falcon\lib\ftb_bc31.lib and this program are linked  */
/* & then make sure that the Options->Compiler->Code generation->Model */
/* is set to Large, and Options->Compiler->Advance Code generation->   */
/* Floating Point is set to 80287, and Options->Compiler->Advance      */
/* Code generation->Instruction Set is set to 80286                    */

/* also I believe, Options->Compiler->Advance Code generation-> Far    */
/* Data Treshold should be set to 65535                                */

/* I found out that c:\falcon\lib\ftb_mscv.lib would create a fatal    */
/* linker error, since it calls for a 32-bit reference not compatible  */
/* with the DOS environment, so don't use it                           */

/* according to Ron Green:                                             */
/* no need to call u_sync and/or u_strip, just rely on default values  */
/* default : u_strip(0)                                                */
/* needs to call u_camera(CAMERA_RS170, GRABULY, FALSE) to prepare for       */
/* external triggering                                                 */
/* external trigger is a high to low transition of a (TTL?) pulse      */

/* Date of insertion of comment: 1998-10-30 19:20:35 US/Indiana time zone
  Question relating to the NUMBER_OF_GROUPS parameter:

Enter number of samples : 1
Voltage range : 25-38

Then, file "output" is created, which contains either ones or zeroes.

If I were to run matlab, and do
  load output
  sum(sum(output))

  what approximate value would the double-dimensioanl sum give?

The reason why I need to know is that I am setting
  MAX_NUMBER_OF_GROUPS and I want it to be on the right order of the
  maximum number of groups that would be created in one sample time
  period.

Answer given: less than one hundred. Thus, setting MAX_NUMBER_OF_GROUPS
to 512 should be plenty.

*/

#define FALSE 0
#define TRUE 1

#define MAX_NUMBER_OF_GROUPS 512

/* Warning: Setting MAX_RECURSION_LEVEL specifies the number of
   nonzero points that may belong to a group.  There may not be
   a group with more than MAX_RECURSION_LEVEL points.
   */
#define MAX_RECURSION_LEVEL 512

#define OUTPUT_FILENAME "output"

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <dos.h>
#include <process.h>
#include <io.h>
#include <time.h>
#include <fcntl.h>
#include <sys\types.h>
#include <sys\timeb.h>
#include <time.h>

#ifndef _MSC_VER
#include <graphics.h>
#include <mem.h>
#include <alloc.h>
#include <sys\stat.h>
#include <dir.h>
#else
#include <bc2msc.h>
#endif

#include <c:\falcon\include\ftbdefs.h>
#include <c:\falcon\include\tms32030.h>
#include <c:\falcon\include\isllib.h>
#include <c:\falcon\include\tmsmsg.h>
#include <c:\falcon\include\image.h>

/*
  Ok, I have found a horrible, awful, ugly bug somewhere in farmalloc.
  if you allocate an array that passes the 65k boundary, then you get
  strange behavior, and everything crashes.  My array of longs is
  4 bytes long, and so you see that 128*128*4 == 2^16, which is the
  boundary I can't seem to pass.  Thus, the limit on DIMX and DIMY
  is 128.  I have put it at 127 just to be safe.  I am looking into
  getting around this, but right now it is a limitation, so don't
  try and change it!! :)

   -- Tim.
*/

#define DIMX 200
#define DIMY 200
#define AREA 40000
//#define AREA DIMX*DIMY
#define XOFFSET 640/2 - DIMX/2
#define YOFFSET 480/2 - DIMY/2

/* INPST bit 0 is VBLANK                                              */
/* INPST bit 1 is FIELD                                               */
/* INPST bit 2 is PLLDET                                              */
/* INPST bit 3 is (unused?)                                           */
/* INPST bit 4 is INP1                                                */
/* INPST bit 5 is INP2                                                */
/* INPST bit 6 is (unused?)                                           */
/* INPST bit 7 is (unused?)                                           */

/* FMASK 0x10 masks off all bits of the Input Status Register except  */
/* bit 4, which is INP1 flag.                                         */
/* in order to connect INP1 to the trigger pulse properly, a wire has */
/* to connect pin 9 of the P2 connector to the incoming trigger pulse */
/* source.                                                            */
/* by default, in Falcon VL board, INP1 is also connected to JP6, and */
/* by default, in Falcon VL board, pin 15 of the P2 connector, which  */
/* receives the external trigger signal, is also connected to JP7.    */
/* thus we use a red wire to jump JP7 to JP6.                         */

/* INP1 is going to be monitored because when it toggles then we know */
/* that a trigger pulse has arrive.                                   */

#define FMASK 0x10

/* we delay the camera triggering by 35 milliseconds with respect to  */
/* the frame grabber trigeering.                                      */
/* I found out (Aug.2, 96) that it takes about 15 milliseconds beyond */
/* the camera triggering to completely transfer the image from the    */
/* Pulnix camera to the UNDERLAY memory of the Falcon VL.             */
/* u_get_pixel should only be called after the said transfer is done. */
/* hence u_get_pixel needs a delay of (35+15=50) milliseconds after   */
/* the frame grabber triggering.                                      */
#define DLY 10

#define ulong unsigned long
#define uint unsigned int
#define uchar unsigned char

void wait_for_picture(void);
int search_matrix (void);
int search_matrix2 (void) ;
int traverse_group (unsigned char x, unsigned char y);
int traverse_group2 (int x, int y);
//int traverse_group_nonrecursive (  unsigned char x, unsigned char y);
int traverse_group_nonrecursive ( int x, int y);

int mass = 0;
int group = 0;

int max_x = DIMX;
int max_y = DIMY;

/*unsigned char probed [DIMX] [DIMY];*/

int group_x [MAX_NUMBER_OF_GROUPS];
int group_y [MAX_NUMBER_OF_GROUPS];
 // These were of type "int" not too long ago..

long weighted_mass [2] = {0};

uchar far *c_buff;
ulong huge *ui_buff;
//uchar far *probed;
uchar probed[AREA] = {0};

FILE *debug;

int main()
{
  FILE *output;
  uchar r[256], g[256], b[256], missed_string[100];
  ulong largest = 0, loc_largest, missed = 0;
  int samples, x, y, lowrange, highrange;
  int i, j;
  uchar videosignal = CAMERA_RS170;

/*  if((probed = (uchar far *) farcalloc(AREA, sizeof (uchar))) == NULL)
  {
    printf("Calloc failed on probed, going away\n");
    exit(1);
  }*/

  if((c_buff = (uchar far *) farmalloc(AREA)) == NULL)
  {
    printf("Malloc failed on c_buff, going away\n");
    exit(1);
  }

  if((ui_buff = (ulong huge *) farcalloc(AREA, sizeof(ulong))) == NULL)
  {
    printf("Calloc failed on ui_buff, going away\n");
    exit(1);
  }

  printf("How many samples do you want to take? ");
  scanf("%d", &samples);
  printf("Enter the low and high range for the voltage: ");
  scanf("%d %d", &lowrange, &highrange);

//  if (( debug = fopen ("debug2","w") ) == NULL) {
//    printf("Unable to open debug file for writing.\n");
//    return 1;
//  };

  if(samples <= 0)
  {
    printf("Invalid entry for number of samples, exiting");
    exit(1);
  }

  if (u_init_falcon(0x69, LoadXilinx_Flag) != OK)
  {
    u_set_accmode(OVERLAY);
    u_set_videomode(3);
    printf("Improper installation or setup of Falcon board\n");
    printf("Unable to continue\n");
    printf("Please enter any alphanumeric value to exit\n");
    getchar();
    exit(1);
  }

  u_camera(videosignal, GRABULY, FALSE);
  u_set_accmode(OVERLAY);

  /* Set the input bias to vary from 0-1.3 V */
  u_inputbias(lowrange, highrange);

  /*
     Clear the OVERLAY memory, to try to get rid
     of that weird garbage we see on the screen
     sometimes
  */
  u_clearframebuf(0, 0, 480);

  /*
     Set the mode back to UNDERLAY, which is where
     we actually want to grab the pictures to
  */
  u_set_accmode(UNDERLAY);

  for(y = 0; y < samples; y++)
  {
    wait_for_picture();

    /* HACK!! */
    //delay(DLY);
    u_get_area(XOFFSET, YOFFSET, DIMX, DIMY, c_buff);

    group = 0;

    search_matrix2 ();

    for (x = AREA; x--;)
      probed [x] = 0; /* for next time */

    for (x = group; x--;)
      {
	ui_buff [group_x [x] * DIMX + group_y [x]] ++;
      }

    //loc_largest = 0;

/*    for(x = 0; x < AREA; x += 10)
    {
      if(c_buff[x] > 0) ui_buff[x] += c_buff[x];
      if(c_buff[x+1] > 0) ui_buff[x+1] += c_buff[x+1];
      if(c_buff[x+2] > 0) ui_buff[x+2] += c_buff[x+2];
      if(c_buff[x+3] > 0) ui_buff[x+3] += c_buff[x+3];
      if(c_buff[x+4] > 0) ui_buff[x+4] += c_buff[x+4];
      if(c_buff[x+5] > 0) ui_buff[x+5] += c_buff[x+5];
      if(c_buff[x+6] > 0) ui_buff[x+6] += c_buff[x+6];
      if(c_buff[x+7] > 0) ui_buff[x+7] += c_buff[x+7];
      if(c_buff[x+8] > 0) ui_buff[x+8] += c_buff[x+8];
      if(c_buff[x+9] > 0) ui_buff[x+9] += c_buff[x+9];
      //if(c_buff[x] > loc_largest) loc_largest = (ulong)c_buff[x];
      //if(ui_buff[x] > largest) largest = ui_buff[x];
    }
*/

    /*
       If all of the pixels in the last picture we grabbed
       were 0, then we must have missed the picture, so
       up the count of missed frames by 1
    */
    //if(loc_largest == 0)
    //{
    //  missed++;
    //}
  }
//  fclose(debug);

  for(x = 0; x < AREA; x++)
  {
    if(ui_buff[x] > largest) largest = ui_buff[x];
  }

  if(largest > 0)
  {
    sprintf(missed_string, "%u", largest);
    u_text(100, 10, 260, missed_string);
    for(x = 0; x < AREA; x++) c_buff[x] = (uchar)(ui_buff[x] / ((double)largest/255.0));
  }

  else
  {
    for(x = 0; x < AREA; x++) c_buff[x] = (uchar)ui_buff[x];
  }

  if ((output = fopen (OUTPUT_FILENAME, "w")) == NULL)
    {
      fprintf (stderr,
	"unable to open output file \"" OUTPUT_FILENAME "\"\n");
      return 1;
    }

  for(x = 0; x < DIMX; x++)
  {
    for(y = 0; y < DIMY; y++)
    {
      fprintf(output, "%u ", ui_buff[x*DIMX + y]);
    }
    fprintf(output, "\n");
  }

  fclose(output);


  /*
     Get the palet from the UNDERLAY, and set to
     OVERLAY mode, and reset the palet in the OVERLAY
     to the palet that was in the UNDERLAY, this give
     us the grayscale palet we were using in the UNDERLAY
     in the OVERLAY as well.

     I disable the UNDERLAY here so we won't be getting
     a signal from the camera anymore.
  */
  u_get_palet(r, g, b);
  u_set_ovmode(DISABLEUNDERLAY);
  u_set_accmode(OVERLAY);
  u_set_palet(r, g, b);

  /*
     Put our composite picture onto the screen
  */
  u_clearframebuf(0, 0, 480);
  u_put_area(0, 0, DIMX, DIMY, c_buff);

  sprintf(missed_string, "Number missed frames: %d", missed);
  u_text(100, 10, 240, missed_string);

  getch();

  u_set_accmode(OVERLAY);
  u_exit_falcon();

  return 0;

}

void wait_for_picture(void)
{
  while(u_get_freg(INPST, FMASK) != 0);
}

int search_matrix (void)
  {
    unsigned char x;
    unsigned char y;

 //   fprintf(debug,"Starting debug file...\n");

//    for (x = 0; x < max_x-1; x+=2) // !!!!!!!!!!!!!!!!!CHANGE BACK!!!!
    for (x = 0; x < max_x-1; x+=2) {  // just did
 //     fprintf(debug,"sync %ld\n",x);
      for (y = max_y; y--;) {
//	fprintf(debug,"considering %d, %d\n",x,y);

// RESTORE THIS -> -> ->	if (c_buff [x * DIMX + y] && !probed [x * DIMX + y])
	if (c_buff [x * DIMX + y])
	  {
//	    fprintf(debug,"[%ld %ld]\n",x,y);
   //	    fflush(debug);
//	    traverse_group_nonrecursive (x, y);
	    traverse_group (x,y);
	    group_x [group] =   (1+(((weighted_mass [0])>>1) / mass))<<1;
	    group_y [group++] = (1+(((weighted_mass [1])>>1) / mass))<<1;
//	    fprintf(debug,"[%ld %ld] -> mass=%ld, group=%ld, group=[%ld %ld]\n",
  //		 x,y,mass,group-1,group_x[group-1],group_y[group-1]);
	    mass = 0;
	    weighted_mass [0] = 0;
	    weighted_mass [1] = 0;
	  }
      }
    }
 //   fflush(debug);
    return 0;
  }

unsigned char ss[MAX_RECURSION_LEVEL]   = {0};
unsigned char xx[MAX_RECURSION_LEVEL] = {0};
unsigned char yy[MAX_RECURSION_LEVEL] = {0};

int traverse_group_nonrecursive ( int x, int y ) {
    int cache, i = 0;
   /* States:
	State     Delta_X   Delta_Y
	 0          -2         0
	 1           0        -1
	 2           2         0
	 3           0         1
    */

  ss[0] = 0; /* begin at state zero */
  xx[0] = x; /* start at (x,y) */
  yy[0] = y;

//  fprintf(debug,"entering main loop..\n");
//  fflush(debug);
//  flushall();
  for(;;) {

 /* {
    struct  time t;

gettime(&t);
//printf("The current time is: %2d:%02d:%02d.%02d\n",
  //     t.ti_hour, t.ti_min, t.ti_sec, t.ti_hund);



   // time_t t;

    // time(&t);
     //printf("Today's date and time: %s\n", ctime(&t));


    fprintf(debug,
       "time{%02d:%02d.%02d}:i=%d, ss[i]=%d, xx[i]=%d, yy[i]=%d, probed=%d, mass=%d\n",
    t.ti_min,t.ti_sec,t.ti_hund, i, ss [i], xx [i], yy [i], probed[xx[i]*DIMX+yy[i]],
	 mass);
    fflush(debug);
    flushall();
    fclose(debug);
    flushall();
    debug=fopen("debug2","a");
  } */
//  if ( i  == 100 ) { return 1;} // woah.

    if (! ss [i]) {
      if (   (xx[i] > max_x-1)
	  || (yy[i] > max_y-1)
	  || probed [xx [i] * DIMX + yy [i]]
	  || ! c_buff [xx [i] * DIMX + yy [i]]  )
      {
	ss[--i]++;
	continue;
      } else {
	mass += (cache=c_buff[xx[i]*DIMX+yy[i]]);
	weighted_mass [0] += xx[i] * cache; weighted_mass [1] += yy[i] * cache;
	probed [xx[i] * DIMX + yy[i]] ++;
	ss [i] ++;
      }
    } else {
      if ( ss [i] == 5 ) {
	if ( !i ) {
	  break;
	}
	ss[--i]++;
	continue;
      } else {
	  ss[i]--;
	  xx[i+1]=((xx[i]+((ss[i]==2)-(!ss[i]))*2));
	  yy[i+1]=((yy[i]+((ss[i]==3)-(ss[i]==1))));
	  ss[i]++; //may be redundant
	  ss[++i] = 0;
      }
    }
  };

  return TRUE;
}

int traverse_group (unsigned char x, unsigned char y)  {
    unsigned char cache;
    if (! (cache = c_buff [x * DIMX + y]) || probed [x * DIMX + y])
      { return FALSE; }
    mass += cache;
    weighted_mass [0] += x * cache; weighted_mass [1] += y * cache;
    probed [x * DIMX + y] ++;

    if (x>=2) traverse_group (x - 2, y);
    if (y) traverse_group (x, y - 1);
    if (x <= max_x - 3) traverse_group (x + 2, y);
    if (y != max_y - 1) traverse_group (x, y + 1);
    return TRUE;
  }


  /* old code:::: */



int search_matrix2 (void)
  {
    int x;
    int y;

 //   fprintf(debug,"Starting debug file...\n");

    for (x = 0; x < max_x-1; x+=2)
    //for (x = max_x; x--;)
      for (y = max_y; y--;)
	if (c_buff [x * DIMX + y] && !probed [x * DIMX + y])
	  {
   //	    fprintf(debug,"ready to enter t_g_nonr\n");
	    traverse_group_nonrecursive (x, y);
	    //traverse_group2 (x,y);
//old	    group_x [group] = (weighted_mass [0] + (mass >> 1)) / mass;
//   old    group_y [group++] = (weighted_mass [1] + (mass >> 1)) / mass;
	    group_x [group] =   (1+(((weighted_mass [0])>>1) / mass))<<1;
	    group_y [group++] = (1+(((weighted_mass [1])>>1) / mass))<<1;


  //   fprintf(debug,"point [%d %d], wm=[%ld %ld], m=%d\n",group_x[group-1],group_y[group-1],weighted_mass[0],weighted_mass[1],mass);
 //    fflush(debug);
	    mass = 0;
	    weighted_mass [0] = 0;
	    weighted_mass [1] = 0;
	  }

    return 0;
  }

int traverse_group2 (int x, int y)
  {
    int cache;

    if (! (cache = c_buff [x * DIMX + y]) || probed [x * DIMX + y])
      {
	return FALSE;
      }
 //   fprintf(debug,"old mass %d, ",mass);
    mass += cache;
//    fprintf(debug,"new mass %d.\n",mass);
    if (x<0  || y<0) {
  0;//    fprintf(debug,"Warning: x,y =[%ld,%ld]\n",x,y);
   };
    weighted_mass [0] += x * cache;
    weighted_mass [1] += y * cache;
    probed [x * DIMX + y] ++;

    if (x>=1)
      traverse_group2 (x - 1, y);
    if (y>=1)
      traverse_group2 (x, y - 1);
    if (x <= max_x - 1)
      traverse_group2 (x + 1, y);
    if (y <= max_y - 1)
      traverse_group2 (x, y + 1);

    return TRUE;
  }
