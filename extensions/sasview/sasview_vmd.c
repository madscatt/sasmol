/*
 * Program to link sasview python methods to C-functions that
 * call imd/vmdsock routines to handle sending and receiving
 * coordinates between sasmol and vmd
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "vmdsock.h"
#include "imd.h"

int send_coordinates_to_vmd(int N,float *x,float *y,float *z,int port,int flag_clear_sock){

  printf(">>> port = %d\n",port);
  fflush(stdout);

  int length,i;
  IMDEnergies energies;
  float coords[3*N];
  float tmp;

  static void *sock;
  static void *clientsock;

  static int Qstart=1;
  if (Qstart)
  {
	  vmdsock_init();
	  sock = vmdsock_create();
	  vmdsock_bind(sock, port);
  	  printf(">>> Waiting for IMD connection on port %d...\n", port);
        vmdsock_listen(sock);
	  clientsock = NULL;
   	  while (!clientsock)
	  {
    		if (vmdsock_selread(sock, 0) > 0)
			{
      			clientsock = vmdsock_accept(sock);
      			if (imd_handshake(clientsock)) clientsock = NULL;
    		}
   	  }
	  sleep(1);
   	  if (vmdsock_selread(clientsock, 0) != 1 || imd_recv_header(clientsock, &length) != IMD_GO) clientsock = NULL;
	  Qstart = 0;
  }

  for ( i = 0 ; i < N ; i++ )
  {
	  coords[3*i] = x[i] ; 
	  coords[3*i+1] = y[i] ; 
	  coords[3*i+2] = z[i] ; 
  }

  sleep(1.0);
  imd_send_energies(clientsock, &energies);
  imd_send_fcoords(clientsock, N, coords); 

  if (flag_clear_sock)
  {
      if (sock) imd_disconnect(clientsock) ;
      vmdsock_destroy(sock); 
  }

  return(0);
}

