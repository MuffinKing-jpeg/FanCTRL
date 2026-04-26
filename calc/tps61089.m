clear("all");

global vin = 6.0:0.25:12;
global Rsw          = 200e3;
global Imax         = 3;
global Rfb_top      = 300e3;
global Rfb_bot      = 33e3;
global Vref         = 1.224;
global Gea          = 190e-6;
global Rsense       = 0.08;
global Cout         = 47e-6;
global Resr         = 2.15e-3;
global L            = 6.8e-6;
global Vout_aprox   = 12;
global n            = 0.85;
global Fc_rounding  = 5000
global Rilim        = 200e3;

Fvout = @(x) (x-Vref)*Rfb_bot/Vref - Rfb_top;

global Vout = fzero(Fvout, Vout_aprox);
global Ro = Vout/Imax;

printf("Input:\n");
printf("Vin min\t=%10.2f V\nVin max\t=%10.2f V\n",min(vin),max(vin));
printf("Rfb top\t=%10.2f kΩ\nRfb bot\t=%10.2f kΩ\n",Rfb_top/1e3,Rfb_bot/1e3);
printf("R fsw\t=%10.2f kΩ\n",Rsw/1e3);
printf("R Ilim\t=%10.2f kΩ\n",Rilim/1e3);
printf("L\t=%10.2f uH\n",L/1e-6);
printf("C out\t=%10.2f uF\n",Cout/1e-6);
printf("R esr\t=%10.2f mΩ\n",Resr/1e-3);
printf("n\t=%10.2f %%\n",n/1e-2);
printf("Imax\t=%10.2f A\n",Imax);

fsw = zeros(size(vin));
fsw0 = 7e5;

Ilim = 1030000/Rilim;

D_arr = zeros(size(vin));

Ipeak_arr = zeros(size(vin));

Frhpz_arr = zeros(size(vin));

function freq = solve_rsw(cur_vin,result0)
  global Rsw;
  Frsw = @(x) 4*(1./x - 86e-9*12.35/cur_vin)/24e-12 - Rsw;
  freq = fzero(Frsw,result0);
end

function D = solve_D(cur_vin)
  global n;
  global Vout;
  D = 1-cur_vin*n/Vout;
end

function Ipeak = solve_Ipeak (vin, vout, fsw)
  global n;
  global Imax;
  global L;
  Idc = vout*Imax/vin*n;
  Ipp = 1/(L*(1/(vout-vin)+1/vin)*fsw);
  Ipeak = Idc + Ipp/2;
end

function Frhpz = solve_Frhpz (D)
  global Ro;
  global L;
  Frhpz = Ro*pow2(1-D)/(2*pi()*L);
end

for k = 1:length(vin)
  fsw(k) = solve_rsw(vin(k),fsw0);
  D_arr(k) = solve_D(vin(k));
  Ipeak_arr(k) = solve_Ipeak(vin(k),Vout, fsw(k));
  Frhpz_arr(k) = solve_Frhpz(D_arr(k));
end

Fc = cast(
    idivide(
      min(min(Frhpz_arr)/5,
          min(fsw)/10),
    cast(Fc_rounding, "uint32"))*Fc_rounding,"double");

R5 = (2*pi()*Vout*Rsense*Fc*Cout)/((1-max(D_arr))*Vref*Gea);
C5 = Ro*Cout/(2*R5);
C6 = Resr*Cout/R5;

printf("\nOutput:\n");
printf("Vout\t=%10.2f V\n", Vout);
printf("Ilim\t=%10.0f A\n", Ilim);
printf("Fc\t=%10.0f Hz\n", Fc);
printf("Ipeak\t=%10.2f A\n", max(Ipeak_arr));
printf("R5\t=%10.0f Ω\n", R5);
printf("C5\t=%10.1e F\n", C5);
printf("C6\t=%10.1e F\n", C6);

printf("\n%7s %11s %9s %11s\n","Vin","Fsw","Ipeak","Duty");
printf("----------------------------------------\n");

for k = 1:length(vin)
  printf("%5.2f V %8.0f Hz %7.2f A %9.4f %%\n",
  vin(k),
  fsw(k),
  Ipeak_arr(k),
  D_arr(k));
end
