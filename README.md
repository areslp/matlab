% تحديد المعلمات الأساسية للترانزستور
VBE = 0.7; % جهد القاعدة-الانبعاث (VBE)
VCE = 5; % جهد الكولكتور-الانبعاث (VCE)
T = [25, 50, 75, 100]; % درجات الحرارة (°C)
beta = 100; % معامل التكبير

% توليد مجموعة من قيم تيارات القاعدة (IB)
IB = linspace(0.1e-3, 10e-3, 100);

% حساب الخواص الستاتيكية للترانزستور لكل درجة حرارة
IC_steady_state = zeros(length(T), length(IB));
for i = 1:length(T)
   IC_steady_state(i, :) = beta * IB;
end

% حساب الخواص الديناميكية للترانزستور لكل درجة حرارة
IC_dynamic = zeros(length(T), length(IB));
for i = 1:length(T)
   IC_dynamic(i, :) = beta * IB + VCE * exp(-VBE / (0.026 * (T(i) + 273.15)));
end

% الرسم البياني للخواص الستاتيكية
figure;
hold on;
for i = 1:length(T)
    plot(IB, IC_steady_state(i, :), 'DisplayName', sprintf('T = %d°C', T(i)));
end
xlabel('تيار القاعدة (IB)');
ylabel('جهد الكولكتور (IC) - خاصية ستاتيكية');
title('تأثير درجة الحرارة على الخواص الستاتيكية للترانزستور الثنائي القطبية');
legend('Location', 'northwest');
grid on;

% الرسم البياني للخواص الديناميكية
figure;
hold on;
for i = 1:length(T)
    plot(IB, IC_dynamic(i, :), 'DisplayName', sprintf('T = %d°C', T(i)));
end
xlabel('تيار القاعدة (IB)');
ylabel('جهد الكولكتور (IC) - خاصية ديناميكية');
title('تأثير درجة الحرارة على الخواص الديناميكية للترانزستور الثنائي القطبية');
legend('Location', 'northwest');
grid on;
matlab
======