# MAKEFILE

UNAME := $(shell uname)

BINDIR = bin
OBJDIR = obj
SRCDIR = src
INCLDIR = src/include

ifeq ($(UNAME), Linux)
# Linux-based
	LDFLAGS = -L/opt/AMDAPP/lib/x86_64 -lOpenCL -lGLEW -lglfw3 -pthread -lGLU -lGL -lrt -lXrandr -lXxf86vm -lXcursor -lXi -lXinerama -lX11 -L/home/sauce/CL_LIB/lib -lMagick++ -lMagickCore
	CCFLAGS = -g -I/opt/AMDAPP/include -I$(INCLDIR) -fopenmp -I/usr/include/ImageMagick -DGLEW_STATIC
endif

ifeq ($(UNAME), Darwin)
# OS-X
	LDFLAGS = -framework OpenCL -framework OpenGL -L/usr/local/homebrew/Cellar/glfw3/3.1.1/lib -L/usr/local/homebrew/Cellar/glew/1.11.0/lib -L/usr/local/homebrew/Cellar/imagemagick/6.9.0-10/lib -lglfw3 -lGLEW -lMagick++-6.Q16 -lMagickCore-6.Q16
	CCFLAGS = -g -framework OpenCL -framework OpenGL -I$(INCLDIR) -I/usr/local/homebrew/Cellar/glew/1.11.0/include -I/usr/local/include -I/usr/local/homebrew/Cellar/imagemagick/6.9.0-10/include/ImageMagick-6
endif

all: $(BINDIR)/test

$(BINDIR)/test: $(OBJDIR)/renderer/main.o \
								$(OBJDIR)/geometry/primitives.o \
								$(OBJDIR)/assets/texture.o \
								$(OBJDIR)/renderer/render_target.o \
								$(OBJDIR)/renderer/camera.o
	g++ -o $@ $^ $(LDFLAGS)

$(OBJDIR)/renderer/main.o: $(SRCDIR)/renderer/main.cpp
	g++ -o $@ -c $< $(CCFLAGS)

$(OBJDIR)/renderer/camera.o: $(SRCDIR)/renderer/camera.cpp
	g++ -o $@ -c $< $(CCFLAGS)

$(OBJDIR)/geometry/primitives.o: $(SRCDIR)/geometry/primitives.cpp
	g++ -o $@ -c $< $(CCFLAGS)

$(OBJDIR)/assets/texture.o: $(SRCDIR)/assets/texture.cpp
	g++ -o $@ -c $< $(CCFLAGS)

$(OBJDIR)/renderer/render_target.o: $(SRCDIR)/renderer/render_target.cpp
	g++ -o $@ -c $< $(CCFLAGS)


.PHONY: all clean

clean:
	rm $(OBJDIR)/renderer/*.o
	rm $(OBJDIR)/assets/*.o
	rm $(OBJDIR)/geometry/*.o
	rm $(BINDIR)/*
