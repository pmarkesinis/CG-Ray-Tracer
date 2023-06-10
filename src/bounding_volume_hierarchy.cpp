#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include "interpolate.h"
#include <glm/glm.hpp>
#include <algorithm>
#include <iostream>

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    Node root;
    float inf = std::numeric_limits<float>::infinity();
    float uppermost = inf;
    float lowermost = -inf;
    root.upperBound = glm::vec3(lowermost);
    root.lowerBound = glm::vec3(uppermost);
    short axis = 0; // 0 - x axis; 1 - y axis; 2 - z axis;
    for (int i = 0; i < m_pScene->meshes.size(); i++) {
        root.children.push_back(-i - 1);
        for (int j = 0; j < m_pScene->meshes[i].triangles.size(); j++) {
            glm::uvec3 triangle = m_pScene->meshes[i].triangles[j];
            root.children.push_back(j);
            updateBounds(m_pScene->meshes[i].vertices[triangle.x], root.lowerBound, root.upperBound);
            updateBounds(m_pScene->meshes[i].vertices[triangle.y], root.lowerBound, root.upperBound);
            updateBounds(m_pScene->meshes[i].vertices[triangle.z], root.lowerBound, root.upperBound);
        }
    }
    if (root.children.size() <= 2) {
        root.isLeaf = true;
        tree.push_back(root);
    } else {
        tree.push_back(root);
        BoundingVolumeHierarchy::constructorHelper(0, axis);
    }
}

void BoundingVolumeHierarchy::constructorHelper(int nodeIdx, int axis)
{
    Node n = tree[nodeIdx];
    if (n.children.size() == 2) {
        n.isLeaf = true;
        glm::uvec3 triangle = m_pScene->meshes[-n.children[0] - 1].triangles[n.children[1]];
        auto v0 = m_pScene->meshes[-n.children[0] - 1].vertices[triangle[0]];
        auto v1 = m_pScene->meshes[-n.children[0] - 1].vertices[triangle[1]];
        auto v2 = m_pScene->meshes[-n.children[0] - 1].vertices[triangle[2]];
        updateBounds(v0, n.lowerBound, n.upperBound);
        updateBounds(v1, n.lowerBound, n.upperBound);
        updateBounds(v2, n.lowerBound, n.upperBound);
        tree[nodeIdx] = n;
        return;
    }
    std::vector<float> values;
    int meshIdx = -1 * n.children[0] - 1; // Currently are saving mesh indexes as negative numbers.
    int triangleIdx;
    int buffer;
    for (int i = 1; i < n.children.size(); i++)
    {
        buffer = n.children[i];
        if (buffer < 0) {
            meshIdx = -buffer - 1;
        } else {
            triangleIdx = buffer;
            glm::uvec3 triangle = m_pScene->meshes[meshIdx].triangles[triangleIdx];
            values.push_back(BoundingVolumeHierarchy::calculateCentroid(m_pScene->meshes[meshIdx].vertices[triangle.x], m_pScene->meshes[meshIdx].vertices[triangle.y], m_pScene->meshes[meshIdx].vertices[triangle.z], axis));
        }
    }
    std::size_t size = values.end() - values.begin();
    std::size_t middleIdx = size / 2;
    std::nth_element(values.begin(), values.begin() + middleIdx, values.end());
    float median = values[middleIdx]; // There are ways to find it without sorting but we could optimize that later
    int medianCount = 0; // How many numbers left to the median are equal to it.
    for (int i = 0; i < values.size() / 2; i++) {
        if (values[i] == median) {
            medianCount++;
        }
    }
    Node left, right;
    left.lowerBound = right.lowerBound = glm::vec3(std::numeric_limits<float>::infinity());
    left.upperBound = right.upperBound = glm::vec3(-std::numeric_limits<float>::infinity());
    meshIdx = -n.children[0] - 1; // Currently are saving mesh indexes as negative numbers.
    triangleIdx = -1;
    bool meshPushedL = false, meshPushedR = false;
    for (int i = 0; i < n.children.size(); i++) {
        buffer = n.children[i];
        if (buffer < 0) {
            meshIdx = -buffer - 1;
            meshPushedL = false;
            meshPushedR = false;
        } else {
            triangleIdx = buffer;
            glm::uvec3 triangle = m_pScene->meshes[meshIdx].triangles[triangleIdx];
            Vertex v0 = m_pScene->meshes[meshIdx].vertices[triangle[0]];
            Vertex v1 = m_pScene->meshes[meshIdx].vertices[triangle[1]];
            Vertex v2 = m_pScene->meshes[meshIdx].vertices[triangle[2]];
            glm::vec3 centroid = BoundingVolumeHierarchy::calculateCentroid(v0, v1, v2);
            if (centroid[axis] < median || medianCount > 0) {
                medianCount--;
                if (!meshPushedL)  {
                    left.children.push_back(-meshIdx - 1);
                    meshPushedL = true;
                }
                left.children.push_back(triangleIdx);
                updateBounds(v0, left.lowerBound, left.upperBound);
                updateBounds(v1, left.lowerBound, left.upperBound);
                updateBounds(v2, left.lowerBound, left.upperBound);
            } else {
                if (!meshPushedR) {
                    right.children.push_back(-meshIdx - 1);
                    meshPushedR = true;
                }   
                right.children.push_back(triangleIdx);
                updateBounds(v0, right.lowerBound, right.upperBound);
                updateBounds(v1, right.lowerBound, right.upperBound);
                updateBounds(v2, right.lowerBound, right.upperBound);
            }
        }
    }
    n.children.clear();
    tree.push_back(left);
    n.children.push_back(tree.size() - 1);
    BoundingVolumeHierarchy::constructorHelper(tree.size() - 1, (axis + 1) % 3);
    tree.push_back(right);
    n.children.push_back(tree.size() - 1);
    BoundingVolumeHierarchy::constructorHelper(tree.size() - 1, (axis + 1) % 3);
    tree[nodeIdx] = n;
}

// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    return helperNumLevels(tree, tree[0]);
}

int BoundingVolumeHierarchy::helperNumLevels(const std::vector<Node> &tree, Node node)
{
    if (node.isLeaf == true)
        return 0;

    return glm::max(helperNumLevels(tree, tree[node.children[0]]), helperNumLevels(tree, tree[node.children[1]])) + 1;
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
int BoundingVolumeHierarchy::numLeaves() const
{
    int count = 0;
    for (Node n : tree) {
        if (n.isLeaf) {
            count++;
        }
    }
    return count;
}

void calculateIndices(std::vector<Node>& nodes, int currentlevel, int idx, int level, std::vector<Node>& tree)
{
    Node node = tree[idx];
    if (level == currentlevel)
        nodes.push_back(node);
    else if (!node.isLeaf) {
        calculateIndices(nodes, currentlevel + 1, node.children[0], level, tree);
        calculateIndices(nodes, currentlevel + 1, node.children[1], level, tree);
    }
}

// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{
    if (level < 0)
        level = 0;

    std::vector<Node> nodes;
    calculateIndices(nodes, 0, 0, level, tree);
    for (int i = 0; i < nodes.size(); i++) {
        Node node = nodes[i];
        AxisAlignedBox box;
        box.lower = node.lowerBound;
        box.upper = node.upperBound;
        drawAABB(box , DrawMode::Wireframe);
    }
}

// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    if (leafIdx > 0 && leafIdx <= tree.size()) {
        int nodeIdx = leafIdx - 1;
        int check = 0; 
        Node triangle;
        for (int i = 0; i < tree.size(); i++) {
            Node node = tree[i];
            if (node.isLeaf && check == nodeIdx) {
                triangle = node;
                break;
            } else if (node.isLeaf) {
                check++;        
            }
        }

        AxisAlignedBox box;
        box.lower = triangle.lowerBound;
        box.upper = triangle.upperBound;

        drawAABB(box, DrawMode::Wireframe);
        Scene* scene = m_pScene;

        int meshInd = -triangle.children[0] - 1;
        int triangleInd = triangle.children[1];
        glm::uvec3 vertex = scene->meshes[meshInd].triangles[triangleInd];
        Vertex v0_ = scene->meshes[meshInd].vertices[vertex.x];
        Vertex v1_ = scene->meshes[meshInd].vertices[vertex.y];
        Vertex v2_ = scene->meshes[meshInd].vertices[vertex.z];
        drawTriangle(v0_, v1_, v2_);
    }
}

void findTriangles(const std::vector<Node>& tree, std::vector<std::pair<int, float>>& set, Ray ray, int idx)
{
    float max = std::numeric_limits<float>::max();
    Node current = tree[idx];
    AxisAlignedBox box;
    box.lower = current.lowerBound;
    box.upper = current.upperBound;
    ray.t = max;


    if (intersectRayWithShape(box, ray)) {
        drawAABB(box, DrawMode::Wireframe, glm::vec3 { 1, 1, 1 });
        float t;
        t = ray.t;
        ray.t = max;
        if (current.isLeaf) {
            set.push_back(std::make_pair(idx, t));
            return;
        }

        AxisAlignedBox box1;
        box1.lower = tree[current.children[0]].lowerBound;
        box1.upper = tree[current.children[0]].upperBound;

        if (intersectRayWithShape(box1, ray)) {
            findTriangles(tree, set, ray, current.children[0]);
        }
        ray.t = max;

        AxisAlignedBox box2;
        box2.lower = tree[current.children[1]].lowerBound;
        box2.upper = tree[current.children[1]].upperBound;

        if (intersectRayWithShape(box2, ray)) {
            findTriangles(tree, set, ray, current.children[1]);
        }
    }
}

float firstIntersection(Ray& ray, const glm::vec3& lower, const glm::vec3& upper) {

    float txmin, tymin, tzmin, txmax, tymax, tzmax;

    txmax = (upper.x - ray.origin.x) / ray.direction.x;
    txmin = (lower.x - ray.origin.x) / ray.direction.x;

    tzmax = (upper.z - ray.origin.z) / ray.direction.z;
    tzmin = (lower.z - ray.origin.z) / ray.direction.z;

    tymax = (upper.y - ray.origin.y) / ray.direction.y;
    tymin = (lower.y - ray.origin.y) / ray.direction.y;

    if (ray.direction.x == 0) {
        txmin = std::numeric_limits<float>::min();
        txmax = std::numeric_limits<float>::max();
    }
    if (ray.direction.y == 0) {
        tymin = std::numeric_limits<float>::min();
        tymax = std::numeric_limits<float>::max();
    }
    if (ray.direction.z == 0) {
        tzmin = std::numeric_limits<float>::min();
        tzmax = std::numeric_limits<float>::max();
    }

    float txout, tyout, tzout;

    txout = glm::max(txmax, txmin);
    tyout = glm::max(tymax, tymin);
    tzout = glm::max(tzmax, tzmin);

    return glm::min(txout, glm::min(tyout, tzout));
}

//void traversalVisualDebug(Ray ray, )




// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h.
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const
{
    // If BVH is not enabled, use the naive implementation.
    if (!features.enableAccelStructure) {
        bool hit = false;
        float t = ray.t;
        Vertex vector0;
        Vertex vector1;
        Vertex vector2;
        // Intersect with all triangles of all meshes.
        for (const auto& mesh : m_pScene->meshes) {
            for (const auto& tri : mesh.triangles) {
                Vertex v0 = mesh.vertices[tri[0]];
                Vertex v1 = mesh.vertices[tri[1]];
                Vertex v2 = mesh.vertices[tri[2]];
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                    if (ray.t < t) {
                        vector0 = v0;
                        vector1 = v1;
                        vector2 = v2;
                        hitInfo.material = mesh.material;
                        t = ray.t;
                    }
                    // if intersects 
                        if (features.enableTextureMapping == true) {
                            // intersection point as input to calculate the barycentric coordinates 
                            glm::vec3 p = ray.origin + ray.direction * ray.t;
                            glm::vec3 barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, p);
                            // calculate the text coords using barycentric coords
                            glm::vec2 texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, barycentricCoord);

                            // update info in hitInfo 
                            hitInfo.texCoord = texCoord;
                            hitInfo.barycentricCoord = barycentricCoord;   
                        }
                    hit = true;
                }
            }
        }

        if (hit) {
            glm::vec3 barycentric = computeBarycentricCoord(vector0.position, vector1.position, vector2.position, ray.origin + ray.direction * ray.t);
            glm::vec3 interpolated = interpolateNormal(vector0.normal, vector1.normal, vector2.normal, barycentric);
            hitInfo.normal = interpolated;
            hitInfo.barycentricCoord = barycentric;
            if (features.enableNormalInterp) {
                drawNormalInterpolation(vector0, vector1, vector2, ray.origin + ray.direction * ray.t, interpolated, ray);
            }

        }
        // Intersect with spheres.
        for (const auto& sphere : m_pScene->spheres)
            hit |= intersectRayWithShape(sphere, ray, hitInfo);
        return hit;
    } else {
        // TODO: implement here the bounding volume hierarchy traversal.
        // Please note that you should use `features.enableNormalInterp` and `features.enableTextureMapping`
        // to isolate the code that is only needed for the normal interpolation and texture mapping features.
        int i, max;
        bool check, hit, draw;
        max = std::numeric_limits<float>::max();
        

        std::vector<std::pair<int, float>> set;
        findTriangles(tree, set, ray, 0);
        std::sort(set.begin(), set.end(), [](const std::pair<int, float>& first, const std::pair<int, float>& last)
            {
            bool res = first.second < last.second;
            return res;
            });
        float finalT = max;
        Vertex v0, v1, v2;
        check = false;
        hit = false;
        draw = true;
        i = 0;
        int j = 1;

        bool result = false;
        for (std::pair<int, float> idx : set) {
            AxisAlignedBox box3;
            box3.lower = tree[idx.first].lowerBound;
            box3.upper = tree[idx.first].upperBound;
            if (finalT < idx.second && hit) {
                if (features.accelerationDataTraversal) {
                    drawAABB(box3, DrawMode::Wireframe, glm::vec3{ 0, 0, 1 });
                }
                continue;
            }
            if (draw == true) {
                drawAABB(box3, DrawMode::Wireframe);
            }

            
            int mesh = -tree[idx.first].children[0] - 1;
            int triangle = tree[idx.first].children[1];

            glm::uvec3& vp = m_pScene->meshes[mesh].triangles[triangle];
            glm::vec3& v0_ = m_pScene->meshes[mesh].vertices[vp.x].position;
            glm::vec3& v1_ = m_pScene->meshes[mesh].vertices[vp.y].position;
            glm::vec3& v2_ = m_pScene->meshes[mesh].vertices[vp.z].position;


            if (intersectRayWithTriangle(v0_, v1_, v2_, ray, hitInfo)) {
                Material& m = m_pScene->meshes[mesh].material;
                hitInfo.normal = glm::cross(v1_ - v0_, v2_ - v0_);
                hitInfo.material = m;
                v0 = m_pScene->meshes[mesh].vertices[vp.x];
                v1 = m_pScene->meshes[mesh].vertices[vp.y];
                v2 = m_pScene->meshes[mesh].vertices[vp.z];
                hit = true;
                check = true;
                i++;
            }
            if (i == j) {
                finalT = firstIntersection(ray, tree[idx.first].lowerBound, tree[idx.first].upperBound);
            }

            if (hit) {
                drawTriangle(v0, v1, v2);
            }

            result = hit;
            
        }
        return result;
    }
}

float BoundingVolumeHierarchy::calculateCentroid(Vertex v1, Vertex v2, Vertex v3, short axis)
{
    if (axis == 0) {
        return (v1.position.x + v2.position.x + v3.position.x) / 3.0f;
    } else if (axis == 1) {
        return (v1.position.y + v2.position.y + v3.position.y) / 3.0f;
    } else if (axis == 2) {
        return (v1.position.z + v2.position.z + v3.position.z) / 3.0f;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

glm::vec3 BoundingVolumeHierarchy::calculateCentroid(Vertex v1, Vertex v2, Vertex v3)
{
    return (v1.position + v2.position + v3.position) / 3.0f;
}

void BoundingVolumeHierarchy::updateBounds(Vertex& v, glm::vec3& lowerBound, glm::vec3& upperBound)
{
    glm::vec3 position = v.position;
    lowerBound.x = std::min(lowerBound.x, position.x);
    lowerBound.y = std::min(lowerBound.y, position.y);
    lowerBound.z = std::min(lowerBound.z, position.z);
    upperBound.x = std::max(upperBound.x, position.x);
    upperBound.y = std::max(upperBound.y, position.y);
    upperBound.z = std::max(upperBound.z, position.z);
}