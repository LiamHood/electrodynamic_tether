function [ S , C ] = StumpfCalc( z ) 
% S is c3 in vallado
% C is c2 in vallado
        if z > 1e-6
            rz = sqrt( z ) ;
            S = ( rz - sin( rz ) )/( sqrt( z^3 ) ) ;
            C = ( 1 - cos( rz ) )/z ;
        elseif z < -1e-6
            rz = sqrt( -z ) ;
            S = ( sinh( rz ) - rz )/( sqrt( (-z)^3 ) ) ;
            C = ( 1 - cosh( rz ) )/( z ) ;
        else 
            S = 1/6 ;
            C = 1/2 ;
        end
end