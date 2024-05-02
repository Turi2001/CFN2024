matricule = 514162;
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
condAtest = cond(A);

x = floor(10*rand(n,1));
b = A*x;

xhat0 = A\b; % solution initiale

nit = 20;
rr = zeros(n,nit+1); eehat = zeros(n,nit+1);
%
% Implémentation d'une boucle de raffinements itératifs for i=1:nit
% et stockage des erreurs et des résidus après chaque raffinement dans rr et ee 
% Donc:
% <initialisations de r, ehat, xhat>
% Initialisation des résidus et des erreurs
rr = zeros(n, nit+1);
eehat = zeros(n, nit+1);

% Initialisation de xhat (solution initiale)
xhat = xhat0;
r =b;
ehat = A\r;
% Boucle de raffinements itératifs
for i = 1:nit
    % Calcul du résidu r = b - A*xhat
    r = r - A*ehat;

    % Résolution du système linéaire A*ehat = r en utilisant la factorisation LU
    [L, U] = lu(A);
    ehat = U \ (L \ r);

    % Mise à jour des résidus et des erreurs
    rr(:, i+1) = r;    % Stocke le résidu r à l'itération i
    eehat(:, i+1) = ehat;    % Stocke l'erreur ehat à l'itération i

    % Mise à jour de xhat pour l'itération suivante
    xhat = xhat + ehat;

    % Note : xhat(i+1) est la nouvelle estimation de la solution à l'itération i+1
end

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
saveas(gcf,"Erreur",'jpeg');