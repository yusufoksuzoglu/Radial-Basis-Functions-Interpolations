% RBF ve Halton Noktaları - Egzersiz Betiği (3D Görselleştirme Eklenmiş Hali)
clear; clc; close all;

% Boyutlar ve nokta katsayıları
dimensions = [1, 2, 3];

for i = 1:length(dimensions)
    d = dimensions(i);
    N = 5^d; % d=1 için 5, d=2 için 25, d=3 için 125 nokta

    fprintf('--- Boyut (d) = %d, Nokta Sayısı (N) = %d ---\n', d, N);

    % 1. Halton Noktalarını Üretme
    p = haltonset(d);
    X = net(p, N); % N adet d-boyutlu Halton noktası
    Z = net(p, 10);
    % 2. Fill-Distance (Dolgu Mesafesi) Hesabı için Değerlendirme Grid'i Oluşturma
    if d == 1
        eval_grid = linspace(0, 1, 100)';
    elseif d == 2
        [x1, x2] = meshgrid(linspace(0, 1, 30));
        eval_grid = [x1(:), x2(:)];
    else
        % 3 boyut için tarama grid'i
        [x1, x2, x3] = meshgrid(linspace(0, 1, 15));
        eval_grid = [x1(:), x2(:), x3(:)];
    end

    % 3. Fill-Distance Hesabı (h_X)
    min_dists = zeros(size(eval_grid, 1), 1);
    for j = 1:size(eval_grid, 1)
        dists = sqrt(sum((X - eval_grid(j,:)).^2, 2));
        min_dists(j) = min(dists);
    end
    fill_distance = max(min_dists);

    fprintf('Hesaplanan Fill-Distance (h_X): %.4f\n\n', fill_distance);

    % 4. Görselleştirme (d = 1, 2 ve 3 için)
    if d == 1
        figure;
        plot(X, zeros(size(X)), 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'r');
        hold on;
        plot(Z, zeros(size(Z)), 'bo', 'MarkerSize', 6, 'LineWidth', 2, 'MarkerFaceColor', 'b');
        xlim([0 1]); ylim([-0.1 0.1]);
        title(sprintf('1D Halton Noktaları'));
        grid on;

    elseif d == 2
        figure;
        plot(X(:,1), X(:,2), 'bo', 'MarkerSize', 6, 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
        xlim([0 1]); ylim([0 1]);
        title(sprintf('2D Halton Noktaları (d=2, N=%d)', N));
        xlabel('x_1'); ylabel('x_2');
        grid on; pbaspect([1 1 1]);

    elseif d == 3
        figure;
        plot3(X(:,1), X(:,2), X(:,3), 'ko', 'MarkerSize', 5, 'LineWidth', 1, 'MarkerFaceColor', 'k');
        xlim([0 1]); ylim([0 1]); zlim([0 1]);
        title(sprintf('3D Halton Noktaları (d=3, N=%d)', N));
        xlabel('x_1'); ylabel('x_2'); zlabel('x_3');
        grid on; box on;
        view(3); % 3D perspektif açısı
    end
end

% RBF Interpolasyonu ve Hata Analizi (Nokta Sayısına Göre Değişim)

n_values = 5:20; % Denenecek nokta sayıları
errors = zeros(size(n_values)); % Hataları saklamak için boş vektör

f_true = @(x) sin(pi * x); % Test fonksiyonu
epsilon = 2.0;             % Şekil parametresi
phi = @(r) exp(-(epsilon * r).^2);
X_test = linspace(0, 1, 200)';
F_real = f_true(X_test);

for idx = 1:length(n_values)
    i = n_values(idx);

    % 1. Halton Noktaları
    p = haltonset(1);
    X = net(p, i); 
    F = f_true(X);

    % 2. RBF Matrisi Kurulumu (A)
    A = zeros(i, i);
    for k = 1:i
        for j = 1:i
            r = norm(X(k) - X(j)); 
            A(k,j) = phi(r);
        end
    end

    % 3. Katsayıları Çözme
    lambda = A \ F;

    % 4. Hata Ölçümü
    S_pred = zeros(size(X_test));
    for k = 1:length(X_test)
        sum_val = 0;
        for j = 1:i
            r = norm(X_test(k) - X(j));
            sum_val = sum_val + lambda(j) * phi(r);
        end
        S_pred(k) = sum_val;
    end

    errors(idx) = max(abs(F_real - S_pred));
    fprintf('Nokta Sayısı (N): %d | Maksimum Hata: %.6f\n', i, errors(idx));
end

% --- ÖNEMLİ ANALİZ GRAFİĞİ ---
% Nokta sayısı arttıkça hatanın nasıl düştüğünü tek bir grafikte görelim:
figure;
semilogy(n_values, errors, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
title('Nokta Sayısına (N) Göre RBF İnterpolasyon Hatası');
xlabel('Nokta Sayısı (N)');
ylabel('Maksimum Hata (Log Ölçek)');
grid on;

Gaussians = zeros(1,1000);
X = linspace(-1,1,1000);
[x1,x2] = meshgrid(linspace(-1,1,1000));
ep = 3;
for i=1:2
    for j=1:2
        if i==1 && j==1
            Gaussians_Laguerre = (3/2 - (ep*X).^ 2).*exp(-(ep*X).^2);
            figure;
            plot(X, Gaussians_Laguerre);
            title('Gaussian-Laguerre (d=1, n=1)');
        end
        if i==1 && j==2
            Gaussians_Laguerre = (15/8 - (5/2).*(ep*X).^2 + (1/2).*(ep*X).^4).*exp(-(ep*X).^2);
            figure;
            plot(X, Gaussians_Laguerre);
            title('Gaussian-Laguerre (d=1, n=2)');
        end
        if i==2 && j==1
            Gaussians_Laguerre = (2 - ((ep*x1).^2 + (ep*x2).^2)).*(exp( - ((ep*x1).^2 + (ep*x2).^2)));
            figure;
            surf(x1,x2,Gaussians_Laguerre);
            shading interp; colormap('jet');            
            title('Gaussian-Laguerre (d=2, n=1)');
        end
        if i==2 && j==2
            Gaussians_Laguerre = (3 - 3*((ep*x1).^2 + (ep*x2).^2) + ...
                0.5*((ep*x1).^2 + (ep*x2).^2).^2) .* ...
                exp(-((ep*x1).^2 + (ep*x2).^2));
            figure;
            surf(x1, x2, Gaussians_Laguerre);
            shading interp; colormap('jet');
            title('Gaussian-Laguerre (d=2, n=2)');
        end
    end
end

%generalized inverse multiquadratics
beta_1 = 1/2;
beta_2 = 1;
ep = 3;
Generalized_inverse = (1 + ((ep*x1).^2 + (ep*x2).^2)).^(-beta_1);
figure;
surf(x1, x2, Generalized_inverse);
shading interp; colormap('jet');
title('Generalized_inverse (d=2, beta = 1/2)');
Generalized_inverse_2 = (1 + ((ep*x1).^2 + (ep*x2).^2)).^(-beta_2);
figure;
surf(x1, x2, Generalized_inverse_2);
shading interp; colormap('jet');
title('Generalized_inverse_2 (d=2, beta = 1)');

Franke = @(x1, x2) 0.75 * exp(-((9*x1 - 2).^2 + (9*x2 - 2).^2) / 4) + ...
    0.75 * exp(-((9*x1 + 1).^2) / 49 - (9*x2 + 1) / 10) + ...
    0.5  * exp(-((9*x1 - 7).^2 + (9*x2 - 3).^2) / 4) - ...
    0.2  * exp(-((9*x1 - 4).^2 + (9*x2 - 7).^2));

eps = 2.0; % Şekil parametresi
phi = @(r) exp(-(eps * r).^2);

i = 50; 

% 1. Halton Noktaları
p = haltonset(2); 
halton_points = net(p, i); 
x1 = halton_points(:,1);
x2 = halton_points(:,2);
F = Franke(x1,x2);

[x_t1,x_t2] = meshgrid(linspace(0, 1, 20));

% 2. RBF Matrisi Kurulumu (A)
A = zeros(i, i);
for k = 1:i
    for j = 1:i
        r = sqrt((x1(k) - x1(j)).^2 + (x2(k) - x2(j)).^2); 
        A(k,j) = phi(r);
    end
end

% 3. Katsayıları Çözme
lambda = A \ F;

% 4. Fonksiyon Elde Edimi
S_pred = zeros(size(x_t1));
for k = 1:numel(x_t1)
    sum_val = 0;
    for j = 1:i
        r = sqrt((x_t1(k) - x1(j)).^2 + (x_t2(k) - x2(j)).^2);
        sum_val = sum_val + lambda(j) * phi(r);
    end
    S_pred(k) = sum_val;
end
F_real = Franke(x_t1, x_t2);
figure;
surf(x_t1, x_t2, F_real);
shading interp; colormap('jet');
title('Franke, d=2');

figure;
surf(x_t1, x_t2, S_pred);
shading interp; colormap('jet');
title('S_pred, d=2');
