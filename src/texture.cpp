#include "texture.h"
#include <framework/image.h>

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    // TODO: implement this function.
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)

   // performing texture lookup
   // resulting color is used in shading to add textures 
  
    // added offset since first pixel centered at (0.5, 0.5)
    int valueI = round(texCoord.x * image.width - 0.5);
    int valueJ = round(texCoord.y * image.height - 0.5);

    // formula: (j * width + i) 
    int index = (image.height - 1 - valueJ) * image.width + valueI;
    return image.pixels[index];
   
}