function params = genParamsFun(Mov, Ref, varargin)
% Set input parser
p = inputParser;
p.CaseSensitive = false;
p.addParameter('mode', 'point2plane');
p.addParameter('Hb_point2plane', logical(1) );  % if this is true, Lie algebra will be used. 
p.addParameter('Mov_Normal', []);
p.addParameter('Ref_Normal', []);
p.addParameter('optimize', 'irls'); 
p.addParameter('Tf0', eye(4)); 
%%%%%%%% P = [1.0 2.0] or [2.0 2.0], both are OK. 
p.addParameter('P', [1.0 2.0]); 
p.addParameter('maxIter_icp', 100);
%%%%% maxIter_em and maxIter_irls are very important. 
p.addParameter('maxIter_em', 100);   % 20
p.addParameter('maxIter_irls', 100); % 10
p.addParameter('tf_eps', 1e-3);
p.addParameter('IS_SHOW', 0);
p.addParameter('verbose', 1);
p.addParameter('bRand', logical(1) );
parser = p;
parser.parse(varargin{:});
params = parser.Results; 
params.Mov = Mov; 
params.Ref = Ref; 
end
