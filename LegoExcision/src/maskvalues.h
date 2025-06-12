/* $Header$ */

/* Values of the mask */

#ifndef MASKVALUES_H
#define MASKVALUES_H

#ifdef CCODE

#define MASK_EXCISED  0.0
#define MASK_BOUNDARY 0.5
#define MASK_ACTIVE   1.0

#endif

#ifdef FCODE

#define MASK_EXCISED  0.0D0
#define MASK_BOUNDARY 0.5D0
#define MASK_ACTIVE   1.0D0

#endif

#endif /* #ifdef MASKVALUES_H */
