#define SPATIAL_DETERMINANT(gxx_,gxy_,gxz_,gyy_,gyz_,gzz_) \
  (-(gxz_)**2*(gyy_) + 2*(gxy_)*(gxz_)*(gyz_) - (gxx_)*(gyz_)**2 - (gxy_)**2*(gzz_) \
   + (gxx_)*(gyy_)*(gzz_))

#define DOTP(gxx_,gxy_,gxz_,gyy_,gyz_,gzz_,x1_,y1_,z1_,x2_,y2_,z2_) \
 ( (gxx_)*(x1_)*(x2_)+(gyy_)*(y1_)*(y2_)+(gzz_)*(z1_)*(z2_)+ \
   (gxy_)*( (x1_)*(y2_)+(y1_)*(x2_) )+(gxz_)*( (x1_)*(z2_)+(z1_)*(x2_) )+\
   (gyz_)*( (y1_)*(z2_)+(z1_)*(y2_) ) )

#define DOTP2(gxx_,gxy_,gxz_,gyy_,gyz_,gzz_,x_,y_,z_)	\
 ( (gxx_)*(x_)**2+(gyy_)*(y_)**2+(gzz_)*(z_)**2+ \
  2.0*( (gxy_)*(x_)*(y_)+(gxz_)*(x_)*(z_)+(gyz_)*(y_)*(z_) ) )
