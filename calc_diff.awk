function abs(m){ 
  return m> 0 ? m : -m;
}
{
    if (NR > 1) {
      print sold, $4, abs($5 - xold),  abs($6 - vold);
    }
    xold = $5;
    vold = $6;
    sold = $3
}

      
