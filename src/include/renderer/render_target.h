#include <stdgl.h>
#include <assets/texture.h>


#ifndef RENDERTARGET_H
#define RENDERTARGET_H

class RenderTarget{
  GLuint fbo;
  BlankTexture* pTex;
  BlankTexture* pDepth;

  int width;
  int height;
  int numSamples;
  bool mipmaps;
public:
  RenderTarget( int width, int height, int, bool );
  RenderTarget( int width, int height, GLenum internalFormat, GLenum format, GLenum type, int, bool );
  void glBind();
  void glBind( GLenum target );
  void glInit();
  void generateMipmaps();
  Texture* getColorTexture();
  Texture* getDepthTexture();
};
#endif
