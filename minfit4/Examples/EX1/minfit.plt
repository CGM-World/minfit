color 2. 1 0 0 
color 1

autolweight 0
lweight 0.5
submargins 0.0 0.0

expand 1.2
tlabel \\oQ1148+384  Z(abs)=0.553361
xlabel \\oRest Frame Velocity [km s\\u-1\\d]
expand 1.0

def pltions
window 1 5 &1
data &2.min
xc 1
yc 2
xlimits
ylimits -0.2 1.6 
ticks 10. 50. 0.1 1.0
expand 0.6
box 0 0 &4 0
expand 1.3
box &3 4 -1 -1 
ltype 1
relocate -1000. 0.
draw 1000. 0.
ltype 0
data &2
xc 2
yc 3
hist
yc 4
hist
color 2 
data &2.min
xc 1
yc 2
connect
color 1
expand 1.0
ylabel \\o &5 &6 &7\\r
expand 1.0
data ticks.dat
yc 1
xc 2
ptype 1 -1 |
points
end

pltions 1 MgI2853   3 -1 Mg I 2853 
pltions 2 MgII2803 -1 -1 Mg II 2803
pltions 3 MgII2796 -1 -1 Mg II 2796
pltions 4 FeII2600 -1 -1 Fe II 2600
pltions 5 FeII2587 -1  5 Fe II 2587

