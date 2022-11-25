/*
 * Stellarium
 * Copyright (C) 2002-2016 Fabien Chereau and Stellarium contributors
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Suite 500, Boston, MA  02110-1335, USA.
 */

/*
  This is the vertex shader for solar system object rendering
 */

ATTRIBUTE highp vec4 vertex; // vertex projected by CPU
ATTRIBUTE mediump vec2 texCoord;
ATTRIBUTE highp vec4 unprojectedVertex; //original vertex coordinate (in km for OBJ models, in AU otherwise)

uniform highp mat4 projectionMatrix;

VARYING mediump vec2 texc; //texture coord
VARYING highp vec3 P; //original unprojected position (in AU)

VARYING highp vec3 normalX;
VARYING highp vec3 normalY;
VARYING highp vec3 normalZ;

void main()
{
    gl_Position = projectionMatrix * vertex;
    texc = texCoord;

    //unprojectedVertex is already in AU
    P = unprojectedVertex.xyz;

    highp vec3 normal = normalize(unprojectedVertex.xyz);
	if (abs(normal.z)==1.0) // avoid invalid pole anomaly
	{
		normalX=vec3(0);
		normalY=vec3(0);
	}
	else
	{
		normalX = normalize(cross(vec3(0,0,1), normal));
		normalY = normalize(cross(normal, normalX));
	}
	normalZ = normal;
}
