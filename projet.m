matricule = 514162
bulletsize = 8; bulletcolor = 'blue';
bullet = {'o','markerfacecolor',bulletcolor,'color',bulletcolor,...
          'markersize',bulletsize};
rand('state',matricule);
if exist('condA') ~= 1, condA = 1.e14; end
if exist('sizeA') ~= 1, sizeA = 10; end
n = sizeA;

% create matrices of size n x n with given condition number
% First make random matrix
X = rand(n,n);
% then find QR-factorisation to obtain an orthogonal matrix U
[U R] = qr(X);
% repeat the procedure to obtain an orthogonal matrix V
X = rand(n,n);
[V R] = qr(X);
% construct a diagonal matrix with condition number as given
Sigma = diag(1+(n-1:-1:0)/(n-1)*(condA-1));
A = U*Sigma*V';
% reduce values of A to the interval [-1,1]
A = A/max(max(abs(A)));
% check that condition number is as given
condAtest = cond(A)

x = floor(10*rand(n,1));
b = A*x;

xhat0 = A\b; % solution initiale

nit = 20;
rr = zeros(nit,1); eehat = zeros(nit,1);
%
% Implémentation d'une boucle de raffinements itératifs for i=1:nit
% et stockage des erreurs et des résidus après chaque raffinement dans rr et ee 
% Donc:
% <initialisations de r, ehat, xhat>
xhat = xhat0;
r = b - A*xhat;
ehat = A\r;  % Calcul de ehat initial

for i = 1:nit
    % Raffinement itératif
    [L, U] = lu(A);  % Factorisation LU de A
    
    % Résoudre le système linéaire pour obtenir la correction Delta x
    delta_x = U \ (L \ r);
    
    % Mise à jour de xhat
    xhat = xhat + delta_x;
    
    % Recalcul du résidu et de ehat
    r = b - A*xhat;
    ehat = A\r;
    
    % Stockage des résidus et erreurs
    rr(i) = norm(r);     % Norme du résidu
    eehat(i) = norm(ehat);  % Norme de l'erreur ehat
end

%for i=1:nit
%    <raffinement itératif: nouvelles valeurs de r, ehat, xhat>
%    rr(i) = r; eehat(i) = ehat;
%end

normrr = sqrt(sum(rr.^2));
normeehat = sqrt(sum(eehat.^2));


figure(1)
normrr = sqrt(sum(rr.^2));
plot((0:nit),log10(normrr),'b','linewidth',3)
hold on
plot((0:nit),log10(normrr),bullet{:})
hold off

figure(2)
normeehat = sqrt(sum(eehat.^2));
plot((0:nit),log10(normeehat),'b','linewidth',3)
hold on
plot((0:nit),log10(normeehat),bullet{:})
hold off