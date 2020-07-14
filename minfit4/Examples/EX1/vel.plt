tlabel q1148c-CWC z=0.553499
define jessica
window 2 7 &1
data &2
xcol 2
ycol 3
xlimits
ylimits -0.1 1.5
box
histo
ycol 4
histo
ylabel &2
ltype 1
relocate -1000. 0.
draw 1000. 0.
ltype 0
end

jessica 2 MgI2853
xlabel Velocity
jessica 1 MgII2796
xlabel Velocity
jessica 3 MgII2803
jessica 5 FeII2600
jessica 7 FeII2587
jessica 9 FeII2383
jessica 11 FeII2374
jessica 13 FeII2344
jessica 4 MnII2606
jessica 6 MnII2594
jessica 8 MnII2577
jessica 10 CaII3970
jessica 12 CaII3935

