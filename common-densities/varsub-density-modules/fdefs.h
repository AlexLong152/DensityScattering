
/** Flag to choose that the code is compile for BGQ. 
 BGQ signals that the code is compile on a BGQ
 memory usage is then measured using the IBM counters  */
#define NOBGQ

/** Flag to choose that the code is compile for standard linux (e.g. JURECA).
 LINUX signals that the code is compile on a linux machine that 
 allows to read the proc file system to get system stats. 
 memory usage is then measured using proc  */
#define NOLINUX

/** Flag to choose that the code is compile for Mac OS.
 MAC signals that the code is compile on a macbook that 
 allows to use getrusage of Debian */

#define MAC


/** Flag chooses debug mode.
 If DEBUG is switched on much more printing is done. Do not used 
 for production runs. */
#define NODEBUG

/** Flag chooses whether unphysical channels are generated.
 If CHUNPHYS is set, Pauli forbidden channels are generated. Can be useful 
 for tests for the code.  */
#define NOCHUNPHYS 


/** in Compton densities, the flag RHODISTR chooses code 
    that stores rho distributed over P12 and uses HDF I/O to write the files 
*/
#define RHODISTR

/** If WFHDF5 is defined, the code uses the HDF5 I/O Routines to read wave functions.
 Otherwise the wave functions in getwave are read in the old formatted way. 
 */
#define WFHDF5

/** If HDF_MPIO is defined, the code uses the MPIO I/O Routines for HDF. 
 Otherwise, the standard driver is used. (usage is defined in transfer 
 properties.) 
 */
#define NOHDF_MPIO

/* need return on last line to allow appending more definitions */

