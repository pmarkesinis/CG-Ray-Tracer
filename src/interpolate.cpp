#include "interpolate.h"
#include <glm/geometric.hpp>

glm::vec3 computeBarycentricCoord (const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{


   glm::vec3 x = v1 - v0;
   glm::vec3 y = v2 - v1;

   float totalArea = glm::length(glm::abs(glm::cross(x, y))) / 2.0f;
   float alpha = glm::length(abs(glm::cross(v2 - v1, p - v1))) / (totalArea * 2.0f);
   float beta = glm::length(abs(glm::cross(v0 - v2, p - v2))) / (totalArea * 2.0f);
   float gamma = glm::length(abs(glm::cross(v1 - v0, p - v0))) / (totalArea * 2.0f);
    

    return glm::vec3 { alpha, beta, gamma };
}

glm::vec3 interpolateNormal(const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 barycentricCoord)
{
    
    float x = n0.x * barycentricCoord.x + n1.x * barycentricCoord.y + n2.x * barycentricCoord.z;
    float y = n0.y * barycentricCoord.x + n1.y * barycentricCoord.y + n2.y * barycentricCoord.z;
    float z = n0.z * barycentricCoord.x + n1.z * barycentricCoord.y + n2.z * barycentricCoord.z;

    glm::vec3 result = glm::vec3{ x, y, z };

    return result;
}

glm::vec2 interpolateTexCoord (const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 barycentricCoord)
{

    glm::vec2 textureCoordinateX = t0 * barycentricCoord.x;
    glm::vec2 textureCoordinateY = t1 * barycentricCoord.y;
    glm::vec2 textureCoordinateZ = t2 * barycentricCoord.z;

    return glm::vec2(textureCoordinateX + textureCoordinateY + textureCoordinateZ);
   
}
