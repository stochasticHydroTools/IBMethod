% ----------------------------------------------------------------------- %
% Title: phiweights
% Author: Jason Kaye, Yuanxun Bao
% Date: June 2013
% Description: Compute one-dimensional kernel weights at Cartesian grid
% points for kernels with up to 6-point support.
%
% Inputs:
% r - Point at which to compute weights.
% kernelID - can be 'stnd3pt', 'stnd4pt', 'bspline3pt', etc.
% K - Parameter for flexible families.
%
% Outputs:
% w - Column vector of weights at Cartesian grid point.
% ----------------------------------------------------------------------- %

function w=phiweights(r,kernelID,K)

w = [0;0;0;0;0;0];

switch kernelID
    case 'bspline4pt'
        
        w(2:5) = bspline4pt([-1-r,-r,1-r,2-r])';
    case 'bspline4pt_d'
        
        w(2:5) = bspline4pt_d([-1-r,-r,1-r,2-r])';
        
    case 'flex6pt'
        
        w = flex6pt([-2-r,-1-r,-r,1-r,2-r,3-r],K)';
        
    case 'flex6pt_d'

        w = flex6pt_d([-2-r,-1-r,-r,1-r,2-r,3-r],K)';
        
    case 'stnd3pt'
        if r <= 1/2
            w(2) = stnd3pt(-1-r);
            w(3) = stnd3pt(-r);
            w(4) = stnd3pt(1-r);
        elseif r>1/2 && r<=1
            w(3) = stnd3pt(-r);
            w(4) = stnd3pt(1-r);
            w(5) = stnd3pt(2-r);
        end
    case 'stnd3pt_d'
        w(2) = stnd3pt_d(-1-r);
        w(3) = stnd3pt_d(-r);
        w(4) = stnd3pt_d(1-r);
        w(5) = stnd3pt_d(2-r);
    case 'new3pt'
        if r <= 1/2
            w(2) = new3pt(-1-r);
            w(3) = new3pt(-r);
            w(4) = new3pt(1-r);
        elseif r>1/2 && r<=1
            w(3) = new3pt(-r);
            w(4) = new3pt(1-r);
            w(5) = new3pt(2-r);
        end
    case 'stnd3ptYang'
        w(2) = stnd3ptYang(-1-r);
        w(3) = stnd3ptYang(-r);
        w(4) = stnd3ptYang(1-r);
        w(5) = stnd3ptYang(2-r);
    case 'stnd4ptYang'
        if r<= 1/2
            w(1) = stnd4ptYang(-2-r);
            w(2) = stnd4ptYang(-1-r);
            w(3) = stnd4ptYang(-r);
            w(4) = stnd4ptYang(1-r);
            w(5) = stnd4ptYang(2-r);
        elseif r>1/2 && r<=1
            w(2) = stnd4ptYang(-1-r);
            w(3) = stnd4ptYang(-r);
            w(4) = stnd4ptYang(1-r);
            w(5) = stnd4ptYang(2-r);
            w(6) = stnd4ptYang(3-r);            
        end
        
    case 'flex5pt'
        
        if r<= 1/2
            w(1) = flex5pt(-2-r,K);
            w(2) = flex5pt(-1-r,K);
            w(3) = flex5pt(-r,K);
            w(4) = flex5pt(1-r,K);
            w(5) = flex5pt(2-r,K);
        elseif r>1/2 && r<=1
            w(2) = flex5pt(-1-r,K);
            w(3) = flex5pt(-r,K);
            w(4) = flex5pt(1-r,K);
            w(5) = flex5pt(2-r,K);
            w(6) = flex5pt(3-r,K);            
        end
    case 'stnd4pt'
        w(2) = stnd4pt(-1-r);
        w(3) = stnd4pt(-r);
        w(4) = stnd4pt(1-r);
        w(5) = stnd4pt(2-r);

    case 'stnd6pt'
        K = 0;
        w(1) = flex6pt(-2-r,K);
        w(2) = flex6pt(-1-r,K);
        w(3) = flex6pt(-r,K);
        w(4) = flex6pt(1-r,K);
        w(5) = flex6pt(2-r,K);
        w(6) = flex6pt(3-r,K);  
        
    case 'new6pt'
        
        K = (59/60)*(1-sqrt(1-(3220/3481)));
        w(1) = flex6pt(-2-r,K);
        w(2) = flex6pt(-1-r,K);
        w(3) = flex6pt(-r,K);
        w(4) = flex6pt(1-r,K);
        w(5) = flex6pt(2-r,K);
        w(6) = flex6pt(3-r,K);

    case 'bspline6pt'

        w = bspline6pt([-2-r,-1-r,-r,1-r,2-r,3-r])';
        
    case 'bspline6pt_d'
        
        w = bspline6pt_d([-2-r,-1-r,-r,1-r,2-r,3-r])';
        
    otherwise
        
        error('Not a valid kernelID.');
        
end
end
