
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


#include <stdio.h>
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
#define DIMX 180
#define DIMY 180
#define AREA DIMX*DIMY
/*
  640/2 - DIMX/2
*/
#define XOFFSET 640/2 - DIMX/2
/*
  480/2 - DIMY/2
*/
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

void main(int argc, char **argv)
{
  FILE *input;
  uchar r[256], g[256], b[256];
  uchar pix_val;
  char line[5000], *token;
  int x = 0, y = 0;

  if(argc != 2)
  {
    printf("Usage: %s <input file>\n", argv[0]);
    exit(1);
  }

  if((input = fopen(argv[1], "r")) == NULL)
  {
    printf("Could not open file %s.\n", argv[1]);
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


  /*
     Get the palet from the UNDERLAY, and set to
     OVERLAY mode, and reset the palet in the OVERLAY
     to the palet that was in the UNDERLAY, this give
     us the grayscale palet we were using in the UNDERLAY
     in the OVERLAY as well.

     I disable the UNDERLAY here so we won't be getting
     a signal from the camera anymore.
  */
//  u_get_palet(r, g, b);
  u_set_ovmode(DISABLEUNDERLAY);
  u_set_accmode(OVERLAY);
  b[0] = r[0] = g[0] = 0;
  for(x = 1; x <= 100; x++)
  {
    b[x] = 255 - x;
    r[x] = x;
    g[x] = x;
  }
  y = 0;
  for(x = 101; x <= 200; x++)
  {
    g[x] = 255 - y;
    r[x] = y;
    b[x] = y;
    y++;
  }
  y = 0;
  for(x = 201; x <= 255; x++)
  {
    r[x] = 255 - y;
    g[x] = y;
    b[x] = y;
    y++;
  }
  u_set_palet(r, g, b);

  x = y = 0;

  while(!feof(input))
  {
    fgets(line, 5000, input);
    line[strlen(line)-1] = '\0';
    if( (token = strtok(line, " ")) != NULL )
    {
      pix_val = (uchar)atoi(token);
      u_put_pixel(x, y, pix_val);
      x++;

      while( (token = strtok(NULL, " ")) != NULL )
      {
	pix_val = (uchar)atoi(token);
	u_put_pixel(x, y, pix_val);
	x++;
      }
    }
    y++;
    x = 0;
  }
  getch();

  u_set_accmode(OVERLAY);
  u_exit_falcon();
}