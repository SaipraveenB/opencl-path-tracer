#include <stdgl.h>
#include <renderer/render_target.h>
#include <iostream>

  RenderTarget::RenderTarget( int width, int height, int numSamples, bool mipmaps ){
    RenderTarget( width, height, GL_RGBA, GL_RGBA, GL_UNSIGNED_BYTE, numSamples, mipmaps );
  }
  RenderTarget::RenderTarget( int width, int height, GLenum internalFormat, GLenum format, GLenum type, int numSamples, bool mipmaps ){
    pTex = new BlankTexture( width, height, internalFormat, format, type, numSamples, mipmaps );
    pDepth = new BlankTexture( width, height, GL_DEPTH_COMPONENT, GL_DEPTH_COMPONENT, GL_FLOAT, numSamples, false );
    this->width = width;
    this->height = height;
    this->numSamples = numSamples;
    this->mipmaps = mipmaps;
    glGenFramebuffers( 1, &fbo );
    glInit();
  }
  Texture* RenderTarget::getColorTexture(){
    return pTex;
  }
  Texture* RenderTarget::getDepthTexture(){
    return pDepth;
  }
  void RenderTarget::generateMipmaps(){
    if( this->mipmaps ){
      pTex->glBind();
      pTex->glGenMipmaps();
    }
  }
  void RenderTarget::glInit(){
    this->glBind();
    pTex->glBind();
    pDepth->glBind();

    if( numSamples ){
      glFramebufferTexture( GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, pTex->glGetInternalTexture() ,0 );
      glFramebufferTexture( GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, pDepth->glGetInternalTexture() ,0 );
    }else{
      glFramebufferTexture2D( GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, pTex->glGetInternalTexture(), 0 );
      glFramebufferTexture2D( GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, pDepth->glGetInternalTexture(), 0 );
    }
  }
  void RenderTarget::glBind(){
    glBindFramebuffer( GL_FRAMEBUFFER, fbo );
    glViewport( 0, 0, width, height );// Set pipeline output to produce data at the right resolution.
  }

  void RenderTarget::glBind( GLenum target ){
    glBindFramebuffer( target, fbo );
    //glViewport( 0, 0, width, height );// Set pipeline output to produce data at the right resolution.
  }
