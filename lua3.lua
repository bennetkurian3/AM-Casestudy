
--Gear parameters 
z1 = ui_numberBox("Number Of Teeth On Gear 1",15);     --Number of tooth in Gear 1
z2 = ui_numberBox("Number Of Teeth On Gear 2 ",18);   --Number of Tooth in Gear 2
m = ui_scalarBox("Module Of Both Gears in mm",6,2);--Gear module
alpha = ui_scalarBox("Pressure Angle of gears in degree ",20,0.1); --pressure angle in degrees                   
h_a_coefficient = ui_scalarBox("Addendum Height Factor",1.0,0.1);   --Addendum Height Factor
h_f_coefficient = ui_scalarBox("Dedendum Height Factor",1.25,0.1); --Dedendum Height Factor
i = ui_scalarBox("Rotation of Gears ",0,10);--rotation factor of both gear
Face_width = ui_numberBox("Face Width of Gear in mm",20) --Gear Face Width
x=ui_scalar("Profile Shift Coefficient",0.5,-1,1) --profile shift factor  

-------------------------------------------
--Assign all math with their Variable
PI= math.pi
sin= math.sin
asin=math.asin
cos= math.cos
tan=math.tan
atan=math.atan
sqrt=math.sqrt
pow=math.pow



--Generate Involute Gear
--Calculating Involute Curve
--r_b is the base radius
--inv_alpha is the involute angle
--inv_c is the involute curve
inv_c = function(r_b, inv_alpha)
inv1= v(r_b *(sin(inv_alpha) - inv_alpha *
cos(inv_alpha)), r_b *(cos(inv_alpha) + inv_alpha*sin(inv_alpha)))
return inv1
end


-----calculating function of angle of involute, rotation, mirror, profile slop

roll_angle_psi = function(psi_b, psi_a) 
R_a=(sqrt((psi_a * psi_a - psi_b * psi_b) / (psi_b * psi_b)))
return R_a
end  

Rotation = function(rotate, coordinate)
R= v(cos(rotate) * coordinate.x + sin(rotate) * coordinate.y,
cos(rotate)*coordinate.y-sin(rotate) *coordinate.x)
return R
end

Mirror = function(coordinate)        
M=v(-coordinate.x, coordinate.y)
return M
end                 

Profile_slop = function(prof)               
P_s= ((prof[5].y - prof[1].y) / (prof[5].x - prof[1].x))
return P_s
end

--------------------------------------------------------- Pressure angle--angle between the two sides of the tooth

function angle_gear(m,x,alpha,z)      
alpha=alpha*PI/180
al=(((PI*m/2) + 2*m*x*tan(alpha))/(z*m/2) + 2*tan(alpha) - 2*alpha) 
return al
end

--calculation of required point on the circle or whole circle

function Circle(a, b, r, angle)
    return v(a + r * cos(angle), b + r * sin(angle))
end 

function extrude(profile, angle_deg, extrusion_direction, scaling_factors, z_steps)
    local n_points = #profile
    local angle_rad = angle_deg / 180 * PI
local vertices = {}
    for j = 0, z_steps - 1 do
 local phi = angle_rad * j / (z_steps -1)
 local vectordirection = extrusion_direction * j / (z_steps - 1)
 local scalefactor = (scaling_factors - v(1, 1, 1)) * (j / (z_steps - 1)) + v(1, 1, 1)
 for i = 1, n_points - 1 do
  vertices[i + j * n_points] = v( vectordirection.x + scalefactor.x * (profile[i].x * cos(phi) - profile[i].y * sin(phi)),vectordirection.y + scalefactor.y * (profile[i].x * sin(phi) + profile[i].y * cos(phi)), vectordirection.z * scalefactor.z)
end
table.insert(vertices, vertices[1 + j * n_points])
end

local vertex_sum_start = v(0, 0, 0)
local vertex_sum_end = v(0, 0, 0)
for i = 1, n_points - 1 
do
vertex_sum_start = vertex_sum_start + vertices[i]
vertex_sum_end = vertex_sum_end + vertices[i + n_points * (z_steps - 1)]
end
table.insert(vertices, vertex_sum_start / (n_points - 1))
    table.insert(vertices, vertex_sum_end / (n_points - 1))

local triangles = {}
local k = 1
for j = 0, z_steps - 2 
do
for i = 0, n_points - 2 
do
triangles[k] = v(i, i + 1, i + n_points) + v(1, 1, 1) * n_points * j
triangles[k + 1] = v(i + 1, i + n_points + 1, i + n_points) + v(1, 1, 1) * n_points * j
k = k + 2
end
end
for i = 0, n_points - 2 do
triangles[k] = v(i + 1, i, n_points * z_steps)
k = k + 1
end
for i = 0, n_points - 2 
do
triangles[k] = v(i + n_points * (z_steps - 1), i + 1 + n_points * (z_steps - 1), n_points * z_steps + 1)
        k = k + 1
end
return polyhedron(vertices, triangles)
end

---Calculation of single tooth in a gear

function gear(z,m,alpha_rad,x,f_r,h_a_coeff,h_f_coeff)
local xy={}  
c_c=0.167
cl = c_c * m                                --clearence of tooth [4]

alpha_rad = alpha * PI / 180                      --pressure angle of gear in radian

D_p = z * m                                            --pitch diameter                                

R_p = D_p / 2                                
 --reference radius

D_base = D_p * cos(alpha_rad)	                       
-- diameter of base circle of the gear 

a= (z1+z2)*m/2 + cl                           
--centre distance      

r_b = D_base / 2                                          
--base radius

d_a = D_p + 2 * m * h_a_coefficient + 2 * m * x              -- addendum diameter (with profile shift) 
 
r_a = d_a / 2                                          --addendum radius

h_a = m * h_a_coefficient                                   --addendum height (h_a_coeff = Addendum Height Factor)

d_f = D_p - 2 * m * h_f_coefficient + 2 * m * x              --dedendum diamter (with profile shift) [1]

r_f = d_f / 2                                           
-- root_radius

h_r = m * h_f_coefficient                                    --dedendum height (h_f_coeff = Dedendum Height Factor)
f_r_c=0.38
f_r = f_r_c * m                                         --Fillet radius [4]

True_dia = sqrt(pow(D_p*sin(alpha_rad)- 2*(h_a-(m*x)-h_r *(1-sin(alpha_rad))),2) + D_base * D_base)	
                                                       --where True dia is true involute diameter
--x = 1 - ( z * math.sin(alpha) * math.sin(alpha) / 2 )     -- profile shift coefficient
True_rad = True_dia / 2                                      --true involute Radius
tooth_angle = m*((PI/2) + 2*x*tan(alpha_rad)) /R_p +2*tan(alpha_rad) -2*alpha_rad
                                                       --Angle between two involute Curve---

Start_SOI = roll_angle_psi(r_b,True_rad)      --Starting of Involute curve
--Start_SOI- Diameter at start of Involute
--End_SOI - Diameter at end of Involute

End_SOI = roll_angle_psi(r_b,r_a)               --ending of involute curve

Points= 15
                                        ----points for Better Accuracy
-------------------------------------

function undercut(m, alpha_rad, h_a_coefficient, h_f_coefficient)
    local h_a = m * h_a_coefficient  -- Addendum height
    local h_f = m * h_f_coefficient  -- Dedendum height

    -- Calculate the undercut
    local undercut = (2 * h_a - h_f) * tan(alpha_rad)

    return undercut
end

-------------------------------------------------

function Fillet_center(prof1, r_f, r_r)
    Slop = (prof1[2].y - prof1[1].y) / (prof1[2].x - prof1[1].x)
 Slop_ang = atan(Slop)

    -- Find the parallel point of involute tangent
 x = prof1[1].x + r_f *cos(Slop_ang + PI / 2)  y = prof1[1].y + r_f * sin(Slop_ang + PI / 2)

 d = (y - Slop * x) / sqrt(Slop * Slop + 1)
     th1 = asin(d / (r_f + r_r)) + Slop_ang
v1=v((r_f + r_r) * cos(th1), (r_f + r_r) * sin(th1))

return v1
end




-------------------------------Calculation For FilletRadius-------------------------------- for Calculation we need Involute Curve to maintain the continuty of curve ---------------------------

local involute = {}                                   
for i = 1, Points do
involute[i] =inv_c(r_b,(Start_SOI + (End_SOI - Start_SOI) * i / Points))
end

Profile_slop_inv = Profile_slop(involute)       pressure_angle = atan(Profile_slop_inv)
center_a = {} 
center_a[1] = Fillet_center(involute,f_r,r_f)
Start_fillet = 2 * PI + atan(center_a[1].y / center_a[1].x)
End_fillet = 3 * PI / 2 + pressure_angle

------------------------starting full gear profile including fillet------------------------

for i=1,z do 
for j=1,Points do                              -- For the Left Fillet
xy[#xy+1]=Rotation(2*PI*i/z,Circle(center_a[1].x,center_a[1].y,f_r,(Start_fillet +(End_fillet-Start_fillet) * j / Points)))
end     

for j=1,Points do                              -- For the Left Involute
xy[#xy+1]=Rotation(2*PI*i/z,inv_c(r_b, (Start_SOI +(End_SOI-Start_SOI) *j / Points)))
end       

for j=Points,1,-1 do                      
-- For the Right involute
xy[#xy+1]=Rotation(2*PI*i/z,Rotation(tooth_angle,Mirror(inv_c(r_b,(Start_SOI +(End_SOI-Start_SOI) *j / Points)))))
end   

for j=Points,1,-1 do                            -- For the Right Fillet
xy[#xy+1]=Rotation(2*PI*i/z,Rotation(tooth_angle,Mirror(Circle(center_a[1].x,center_a[1].y,f_r,(Start_fillet +(End_fillet-Start_fillet) * j / Points)))))
end 

end
xy[#xy+1]=xy[1]  
return xy                               
end

--Gear 1
Gear1 = gear(z1,m,alpha_rad,x,f_r,h_a_coeff,h_f_coeff)          
--assigning the gear tooth values to Gear1
gear_angle= angle_gear(m,x,alpha,z1)            --assigning the angle_gear function to gear_angle which will allow to calculate the angle between the two sides of the tooth
rotation_1=rotate(0,0,gear_angle*90/PI)
rotation_2=rotate(0,0,i)                 
--to rotate  gears in z axis with rotation_2
Evoloid_gear = extrude(Gear1 ,0, v(0,0,-Face_width), v(1,1,1), 20) 
--extruding the Gear1 
shaft1 = ccylinder(5,Face_width+20)           -- shaft hole of Gear1
Gear1 = difference(Evoloid_gear,shaft1)
-------Involute Gear1 assembly pins, shafts, and handle-------------------------

base = cylinder(7.5,5)
h1 = cylinder(1.4,3)
base = difference(base,h1)  --base for the pin
p1 = cylinder(4.9,25)
p2 = cylinder(2.5,25)
p = difference(p1,p2)
pin=union(base,p)
emit(translate(0,0,-25)*pin)
emit(translate(0,0,0)*rotation_2*rotation_1*Gear1,7)          
--creating the Gear1 which will then mate with the Gear2 

--Gear 2 (Yellow colour)-----------------------------------------------------------------

Gear2 = gear(z2,m,alpha_rad,-x,f_r,h_a_coeff,h_f_coeff)          --assigning the gear tooth values to Gear

gear_angle=angle_gear(m,-x,alpha,z2)                             

rotation_1=rotate(0,0,-180 -(360/z2-gear_angle*180/PI)/2)
rotation_2=rotate(0,0,-i*z1/z2)                                  --to rotate  gears in z axis with rotation_2

Gear2 = extrude(Gear2, 0, v(0, 0,-Face_width), v(1, 1, 1), 20)  -- Extrude Gear2 

shaft2 = ccylinder(5, Face_width + 20)  -- Shaft hole of Gear2

Gear2 = difference(Gear2, shaft2)

-- Involute Gear2 assembly pins, shafts, and handle
base = cylinder(11, 5)
h1 = cylinder(1.4, 3)
base = difference(base, h1)  -- Base for pin
p1 = cylinder(4.9, 25)
p2 = cylinder(2.4, 25)
p = difference(p1, p2)
pin = union(base, p)
emit(translate(0, a, -25) * pin)

-- Creating a font object for text
f = font()
text = f:str('333', 1.5)

-- Printing the text on Gear2 (commented out)
-- emit(translate(-20, a + 15, 0) * scale(2, 2, 1.5) * text, 4)

-- Shaft for the handle
c1 = translate(25, a, 0) * cylinder(5, 25)

-- Rounding of shaft with sphere
s1 = translate(25, a, 25) * sphere(5)

Handle = translate(0, 0, 0) * union(c1, s1)  -- Handle on Gear2

-- Creating Gear2 which will mate with Gear1
emit(translate(0, a, 0) * rotation_2 * rotation_1 * Gear2, 5)
