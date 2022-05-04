function Gfm_ctrlstate = estimate_lagrange(t_vec,Gf_vec)

tk      = t_vec(3);
tkm1    = t_vec(2);
tkm2    = t_vec(1);
Gf      = Gf_vec(3);
Gfm1    = Gf_vec(2);
Gfm2    = Gf_vec(1);

Gfm_ctrlstate = ( (tk - tkm1)/( (tkm2-tkm1)*(tkm2-tk) ) ) * Gfm2 ...
    + ( (tk - tkm2)/( (tkm1-tkm2)*(tkm1-tk) ) ) * Gfm1 ...
    + ( (2*tk - tkm2 - tkm1)/( (tk-tkm1)*(tk-tkm2) ) ) * Gf;

end