# MAKEFILE

UNAME := $(shell uname)

BINDIR = bin
OBJDIR = obj
SRCDIR = src
INCLDIR = src/include

ifeq ($(UNAME), Linux)
# Linux-based
	LDFLAGS = -L/opt/AMDAPP/lib/x86_64 -lOpenCL
	CCFLAGS = -I$(INCLDIR)
endif

ifeq ($(UNAME), Darwin)
# OS-X
	LDFLAGS = -framework OpenCL -framework OpenGL -L/usr/local/homebrew/Cellar/glfw3/3.1.1/lib -L/usr/local/homebrew/Cellar/glew/1.11.0/lib -L/usr/local/homebrew/Cellar/imagemagick/6.9.0-10/lib -lglfw3 -lGLEW -lMagick++-6.Q16 -lMagickCore-6.Q16
	CCFLAGS = -framework OpenCL -framework OpenGL -I$(INCLDIR) -I/usr/local/homebrew/Cellar/glew/1.11.0/include -I/usr/local/include -I/usr/local/homebrew/Cellar/imagemagick/6.9.0-10/include/ImageMagick-6
endif

all: $(BINDIR)/test

$(BINDIR)/test: $(OBJDIR)/renderer/main.o \
								$(OBJDIR)/geometry/object.o \
								$(OBJDIR)/assets/texture.o \
								$(OBJDIR)/renderer/render_target.o
	g++ -o $@ $^ $(LDFLAGS)

$(OBJDIR)/renderer/main.o: $(SRCDIR)/renderer/main.cpp

	g++ -o $@ -c $< $(CCFLAGS)

$(OBJDIR)/geometry/object.o: $(SRCDIR)/geometry/object.cpp

	g++ -o $@ -c $< $(CCFLAGS)

$(OBJDIR)/assets/texture.o: $(SRCDIR)/assets/texture.cpp

	g++ -o $@ -c $< $(CCFLAGS)

$(OBJDIR)/renderer/render_target.o: $(SRCDIR)/renderer/render_target.cpp

	g++ -o $@ -c $< $(CCFLAGS)


.PHONY: all clean

clean:
	rm *.o
	rm test
