% t, sec	V, m/s	z, km	rho, kg m-3	  p, mb		T, K    -dz/dt, m/s
clear
vdat = textscan(fopen('viking.txt'),'%d%f%f%f%f%f%f','headerlines',1);


zs{1} = vdat{3};
rho{1} = vdat{4};
pres{1} = vdat{5};
temp{1} = vdat{6};

vdat2 = textscan(fopen('viking2.txt'),'%f%f%f%f%f%f','headerlines',1);


zs{end+1} = vdat2{3};
rho{end+1} = vdat2{4};
pres{end+1} = vdat2{5};
temp{end+1} = vdat2{6};

close all

for i = 1:2
    figure(1)
    hold on
    scatter(zs{i},log(rho{i}));
    title("Density vs Height");
    ylabel("log(Density in [kg/m^3])")
    xlabel("Height above Martian surface [km]")

    figure(2)
    hold on
    scatter(zs{i},log(pres{i}));
    title("Pressure vs Height");
    ylabel("log(Pressure in millibar)")
    xlabel("Height above Martian surface [km]")
    
    figure(3)
    hold on
    scatter(zs{i},temp{i});
    title("Temperature vs Height");
    ylabel("Temperature in Kelvin")
    xlabel("Height above Martian surface [km]")
end