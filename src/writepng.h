#ifndef WRITE_PNG_H
#define WRITE_PNG_H

// write a png file named file_name of size width x height
// data has to be a 4*width*height unsigned char (byte) array of
// RGBA data.
int writePNG(const char *file_name, int width, int height, unsigned char *data);

#endif
