from=$SCARF:/work4/projects/sciml/sciml_projects/kdl/hybrid_ulvz/simulations/outputs
file=output/stations/ARRAY/axisem3d_synthetics.nc
to=outputs/nc/$1
mkdir -p $to
rsync -rP $from/$1/04_solve_prem/$file $to/prem.nc
rsync -rP $from/$1/06_solve_1d/$file $to/1d.nc
rsync -rP $from/$1/08_solve_2d/$file $to/2d.nc
rsync -rP $from/$1/12_solve_3d/$file $to/3d.nc
