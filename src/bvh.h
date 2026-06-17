/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Alexander Hanke
--------------------------------------------------------------------*/

#ifndef BOUNDARY_VOLUME_HIERARCHY_H_
#define BOUNDARY_VOLUME_HIERARCHY_H_

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <cmath>
#include <limits>

struct geo_triangle
{
public:
    geo_triangle(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3, const Eigen::Vector3d& n) : p1(p1), p2(p2), p3(p3), n(n) {};

    Eigen::Vector3d p1;
    Eigen::Vector3d p2;
    Eigen::Vector3d p3;

    Eigen::Vector3d n;
};

// Primitive types: triangle, quadrilateral

class primitive
{
public:
    virtual ~primitive() = default;
    virtual void add_triangles(std::vector<geo_triangle>& triangles) {};
    std::pair<Eigen::Vector3d, Eigen::Vector3d> get_bounding_box() const {return std::make_pair(min, max);};

protected:
    primitive() = default;
    primitive(const Eigen::Vector3d& min, const Eigen::Vector3d& max) : min(min), max(max) {}

    Eigen::Vector3d min;
    Eigen::Vector3d max;
};


class box final : public primitive
{
public:
    box() = delete;
    box(const Eigen::Vector3d& min, const Eigen::Vector3d& max) : primitive(min, max) {}
    ~box() = default;

    void add_triangles(std::vector<geo_triangle>& triangles)
    {
        triangles.emplace_back(Eigen::Vector3d(max.x(), min.y(), min.z()), Eigen::Vector3d(max.x(), max.y(), min.z()), Eigen::Vector3d(max.x(), max.y(), max.z()), Eigen::Vector3d(1.0, 0.0, 0.0));
        triangles.emplace_back(Eigen::Vector3d(max.x(), min.y(), min.z()), Eigen::Vector3d(max.x(), max.y(), max.z()), Eigen::Vector3d(max.x(), min.y(), max.z()), Eigen::Vector3d(1.0, 0.0, 0.0));

        triangles.emplace_back(Eigen::Vector3d(max.x(), max.y(), min.z()), Eigen::Vector3d(min.x(), max.y(), min.z()), Eigen::Vector3d(min.x(), max.y(), max.z()), Eigen::Vector3d(0.0, 1.0, 0.0));
        triangles.emplace_back(Eigen::Vector3d(max.x(), max.y(), min.z()), Eigen::Vector3d(min.x(), max.y(), max.z()), Eigen::Vector3d(max.x(), max.y(), max.z()), Eigen::Vector3d(0.0, 1.0, 0.0));

        triangles.emplace_back(Eigen::Vector3d(min.x(), min.y(), max.z()), Eigen::Vector3d(max.x(), min.y(), max.z()), Eigen::Vector3d(max.x(), max.y(), max.z()), Eigen::Vector3d(0.0, 0.0, 1.0));
        triangles.emplace_back(Eigen::Vector3d(min.x(), min.y(), max.z()), Eigen::Vector3d(max.x(), max.y(), max.z()), Eigen::Vector3d(min.x(), max.y(), max.z()), Eigen::Vector3d(0.0, 0.0, 1.0));

        triangles.emplace_back(Eigen::Vector3d(min.x(), max.y(), min.z()), Eigen::Vector3d(min.x(), min.y(), min.z()), Eigen::Vector3d(min.x(), min.y(), max.z()), Eigen::Vector3d(-1.0, 0.0, 0.0));
        triangles.emplace_back(Eigen::Vector3d(min.x(), max.y(), min.z()), Eigen::Vector3d(min.x(), min.y(), max.z()), Eigen::Vector3d(min.x(), max.y(), max.z()), Eigen::Vector3d(-1.0, 0.0, 0.0));

        triangles.emplace_back(Eigen::Vector3d(min.x(), min.y(), min.z()), Eigen::Vector3d(max.x(), min.y(), min.z()), Eigen::Vector3d(max.x(), min.y(), max.z()), Eigen::Vector3d(0.0, -1.0, 0.0));
        triangles.emplace_back(Eigen::Vector3d(min.x(), min.y(), min.z()), Eigen::Vector3d(max.x(), min.y(), max.z()), Eigen::Vector3d(min.x(), min.y(), max.z()), Eigen::Vector3d(0.0, -1.0, 0.0));

        triangles.emplace_back(Eigen::Vector3d(min.x(), min.y(), min.z()), Eigen::Vector3d(max.x(), min.y(), min.z()), Eigen::Vector3d(max.x(), max.y(), min.z()), Eigen::Vector3d(0.0, 0.0, -1.0));
        triangles.emplace_back(Eigen::Vector3d(min.x(), min.y(), min.z()), Eigen::Vector3d(max.x(), max.y(), min.z()), Eigen::Vector3d(min.x(), max.y(), min.z()), Eigen::Vector3d(0.0, 0.0, -1.0));
    };
};

class cylinder final : public primitive
{
public:
    cylinder() = delete;
    cylinder(const Eigen::Vector3d& base_center, const Eigen::Vector3d& axis, double radius, double height)
        : base_center(base_center), axis(axis.normalized()), radius(radius), height(height)
    {
        // Conservative AABB: each end-cap circle of radius r with normal a_hat
        // extends r*sqrt(1 - a_i^2) along dimension i.
        Eigen::Vector3d a = this->axis;
        Eigen::Vector3d r_ext(radius * std::sqrt(1.0 - a.x()*a.x()),
                               radius * std::sqrt(1.0 - a.y()*a.y()),
                               radius * std::sqrt(1.0 - a.z()*a.z()));
        Eigen::Vector3d C_top = base_center + height * a;
        min = base_center.cwiseMin(C_top) - r_ext;
        max = base_center.cwiseMax(C_top) + r_ext;
    }
    ~cylinder() = default;

    void add_triangles(std::vector<geo_triangle>& triangles)
    {
        const int N = 64;
        const double dPhi = 2.0 * M_PI / N;

        // Orthonormal frame perpendicular to axis
        Eigen::Vector3d a = axis; // already normalized in ctor
        Eigen::Vector3d ref = (std::abs(a.x()) < 0.9) ? Eigen::Vector3d(1, 0, 0) : Eigen::Vector3d(0, 1, 0);
        Eigen::Vector3d u = a.cross(ref).normalized();
        Eigen::Vector3d v = a.cross(u).normalized();

        Eigen::Vector3d C_base = base_center;
        Eigen::Vector3d C_top  = base_center + height * a;

        for (int n = 0; n < N; ++n)
        {
            double phi  = n * dPhi;
            double phi2 = phi + dPhi;

            Eigen::Vector3d r0 = radius * (std::cos(phi)  * u + std::sin(phi)  * v);
            Eigen::Vector3d r1 = radius * (std::cos(phi2) * u + std::sin(phi2) * v);

            // Bottom cap: reversed winding gives outward normal = -a
            triangles.emplace_back(C_base, C_base + r1, C_base + r0, -a);

            // Top cap: forward winding gives outward normal = +a
            triangles.emplace_back(C_top, C_top + r0, C_top + r1, a);

            // Side: average radial direction at the segment mid-angle as face normal
            Eigen::Vector3d n_side = (r0 + r1).normalized();

            // 1st side triangle: base@phi, top@phi2, base@phi2
            triangles.emplace_back(C_base + r0, C_top + r1, C_base + r1, n_side);

            // 2nd side triangle: base@phi, top@phi, top@phi2
            triangles.emplace_back(C_base + r0, C_top + r0, C_top + r1, n_side);
        }
    }

private:
    Eigen::Vector3d base_center;
    Eigen::Vector3d axis;
    double radius;
    double height;
};

class wedge final : public primitive
{
public:
    wedge() = delete;
    wedge(const Eigen::Vector3d& base_center, const Eigen::Vector3d& axis, double width, double height) : base_center(base_center), axis(axis), width(width), height(height) {}
    ~wedge() = default;

private:
    Eigen::Vector3d base_center;
    Eigen::Vector3d axis;
    double width;
    double height;
};

class hexahedron final : public primitive
{public:
    hexahedron() = delete;
    hexahedron(const Eigen::Vector3d& base_center, const Eigen::Vector3d& axis, double width, double height, double depth) : base_center(base_center), axis(axis), width(width), height(height), depth(depth) {}
    ~hexahedron() = default;

private:
    Eigen::Vector3d base_center;
    Eigen::Vector3d axis;
    double width;
    double height;
    double depth;
};

class stl_object final : public primitive
{
public:
    stl_object() = delete;
    stl_object(const std::vector<geo_triangle>& triangles) : triangles(triangles)
    {
        min = Eigen::Vector3d::Constant( std::numeric_limits<double>::max());
        max = Eigen::Vector3d::Constant(-std::numeric_limits<double>::max());
        for (const auto& t : triangles)
        {
            min = min.cwiseMin(t.p1).cwiseMin(t.p2).cwiseMin(t.p3);
            max = max.cwiseMax(t.p1).cwiseMax(t.p2).cwiseMax(t.p3);
        }
    }
    ~stl_object() = default;

private:
    std::vector<geo_triangle> triangles;
};

class bvh_node
{
public:
    bvh_node() = default;
    ~bvh_node() = default;

private:
    box bounding_box;
    std::vector<primitive*> primitives;
    std::vector<bvh_node*> children;
};

class bvh
{

};

#endif