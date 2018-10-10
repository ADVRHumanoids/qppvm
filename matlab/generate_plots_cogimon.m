close all
clc


% left_leg = 7:12;
% right_leg = left_leg + 6;

legs = {left_leg, right_leg, left_arm, right_arm};
leg_names = {'Left leg', 'Right leg','Left arm', 'Right arm'};

[b,a] = butter(2,0.01,'low');
tau_m_filt = filtfilt(b,a,joint_effort')';

xlims = [0 size(joint_effort, 2)/100];

time = (1:size(joint_effort, 2))/100;

for i = 1 : 2

    fig = figure;
    fig.Units = 'centimeters';
    fig.Position = [10 10 25 17];
    fig.Color = [1 1 1];
    
    subplot(2,2,1)
%     plot(time, joint_effort(legs{i}(1), :)')
    xlim(xlims)
    xlabel('Time [s]')
    ylabel('Torque [Nm]')
    title([leg_names{i}, ' Hip roll'])
    hold on
    grid on
    plot(time, effort_reference(legs{i}(1), :)' - tau_m_filt(legs{i}(1), :)', 'LineWidth', 3)
%     plot(time, effort_reference(legs{i}(1), :)')
%     plot(time, gcomp_raw(legs{i}(1), :)')
%     plot(time, gcomp_imu(legs{i}(1), :)') 
    
    subplot(2,2,2)
%     plot(time, joint_effort(legs{i}(2), :)')
    xlim(xlims)
    xlabel('Time [s]')
    ylabel('Torque [Nm]')
    title([leg_names{i}, ' Hip pitch'])
    hold on
    grid on
    plot(time, effort_reference(legs{i}(2), :)' - tau_m_filt(legs{i}(2), :)', 'LineWidth', 3)
%     plot(time, effort_reference(legs{i}(2), :)')
%     plot(time, gcomp_raw(legs{i}(2), :)')
%     plot(time, gcomp_imu(legs{i}(2), :)')
    
    subplot(2,2,3)
%     plot(time, joint_effort(legs{i}(3), :)')
    xlim(xlims)
    xlabel('Time [s]')
    ylabel('Torque [Nm]')
    title([leg_names{i}, ' Hip yaw'])
    hold on
    grid on
    plot(time, effort_reference(legs{i}(3), :)' - tau_m_filt(legs{i}(3), :)', 'LineWidth', 3)
%     plot(time, effort_reference(legs{i}(3), :)')
%     plot(time, gcomp_raw(legs{i}(3), :)')
%     plot(time, gcomp_imu(legs{i}(3), :)')
    
    subplot(2,2,4)
%     plot(time, joint_effort(legs{i}(4), :)')
    xlim(xlims)
    xlabel('Time [s]')
    ylabel('Torque [Nm]')
    title([leg_names{i}, ' Knee pitch'])
    hold on
    grid on
    plot(time, effort_reference(legs{i}(4), :)' - tau_m_filt(legs{i}(4), :)', 'LineWidth', 3)
%     plot(time, effort_reference(legs{i}(4), :)')
%     plot(time, gcomp_raw(legs{i}(4), :)')
%     plot(time, gcomp_imu(legs{i}(4), :)')
    legend({'Measured torque', 'Measured torque (filtered)', 'Gcomp', 'Cgomp (IMU compensated)'}, 'Location', 'best')
    
%     fig = figure;
%     fig.Units = 'centimeters';
%     fig.Position = [10 10 25 17];
%     fig.Color = [1 1 1];
%     
%     plot(time, tau_m_filt(legs{i}(1:4),:)' - gcomp_imu(legs{i}(1:4),:)', 'LineWidth', 2)
%     xlim(xlims)
%     xlabel('Time [s]')
%     ylabel('Torque [Nm]')
%     title([leg_names{i}, ' Torque offset']);
%     hold on
%     grid on
%     legend({' Hip pitch', ' Hip pitch', ' Hip yaw', ' Knee pitch'}, 'Location', 'best')
    
    savefig(['/tmp/torque_offsets_'  leg_names{i} '.fig']);
    print(['/tmp/torque_offsets_'  leg_names{i} '.eps'], '-depsc');
    

end

for i = 3 : 4

    fig = figure;
    fig.Units = 'centimeters';
    fig.Position = [10 10 25 17];
    fig.Color = [1 1 1];
    
    subplot(3,3,1)
%     plot(time, effort_reference(legs{i}(1), :)'-joint_effort(legs{i}(1), :)')
    xlim(xlims)
    xlabel('Time [s]')
    ylabel('Torque [Nm]')
    title([leg_names{i}, ' 1'])
    hold on
    grid on
    plot(time, effort_reference(legs{i}(1), :)'- tau_m_filt(legs{i}(1), :)', 'LineWidth', 3)
%     plot(time, effort_reference(legs{i}(1), :)')
%     plot(time, gcomp_raw(legs{i}(1), :)')
%     plot(time, gcomp_imu(legs{i}(1), :)') 
    
    subplot(3,3,2)
%     plot(time,  joint_effort(legs{i}(2), :)')
    xlim(xlims)
    xlabel('Time [s]')
    ylabel('Torque [Nm]')
    title([leg_names{i}, ' 2'])
    hold on
    grid on
    plot(time, effort_reference(legs{i}(2), :)'- tau_m_filt(legs{i}(2), :)', 'LineWidth', 3)
%     plot(time, effort_reference(legs{i}(2), :)')
%     plot(time, gcomp_raw(legs{i}(2), :)')
%     plot(time, gcomp_imu(legs{i}(2), :)')
    
    subplot(3,3,3)
%     plot(time, effort_reference(legs{i}(3), :)'-joint_effort(legs{i}(3), :)')
    xlim(xlims)
    xlabel('Time [s]')
    ylabel('Torque [Nm]')
    title([leg_names{i}, ' 3'])
    hold on
    grid on
     plot(time, effort_reference(legs{i}(3), :)'- tau_m_filt(legs{i}(3), :)', 'LineWidth', 3)
%     plot(time, effort_reference(legs{i}(3), :)')
%     plot(time, gcomp_raw(legs{i}(3), :)')
%     plot(time, gcomp_imu(legs{i}(3), :)')
    
    subplot(3,3,4)
%     plot(time, joint_effort(legs{i}(4), :)')
    xlim(xlims)
    xlabel('Time [s]')
    ylabel('Torque [Nm]')
    title([leg_names{i}, ' 4'])
    hold on
    grid on
    plot(time, effort_reference(legs{i}(4), :)'-tau_m_filt(legs{i}(4), :)', 'LineWidth', 3)
%     plot(time, effort_reference(legs{i}(4), :)')
%     plot(time, gcomp_raw(legs{i}(4), :)')
%     plot(time, gcomp_imu(legs{i}(4), :)')
    legend({'Measured torque', 'Measured torque (filtered)', 'Gcomp', 'Cgomp (IMU compensated)'}, 'Location', 'best')
    
    subplot(3,3,5)
%     plot(time, joint_effort(legs{i}(5), :)')
    xlim(xlims)
    xlabel('Time [s]')
    ylabel('Torque [Nm]')
    title([leg_names{i}, ' 5'])
    hold on
    grid on
    plot(time, effort_reference(legs{i}(5), :)'- tau_m_filt(legs{i}(5), :)', 'LineWidth', 3)
%     plot(time, effort_reference(legs{i}(5), :)')
%     plot(time, gcomp_raw(legs{i}(1), :)')
%     plot(time, gcomp_imu(legs{i}(1), :)') 
    
    subplot(3,3,6)
%     plot(time, joint_effort(legs{i}(6), :)')
    xlim(xlims)
    xlabel('Time [s]')
    ylabel('Torque [Nm]')
    title([leg_names{i}, ' 6'])
    hold on
    grid on
    plot(time, effort_reference(legs{i}(6), :)'-tau_m_filt(legs{i}(6), :)', 'LineWidth', 3)
%     plot(time, effort_reference(legs{i}(6), :)')
%     plot(time, gcomp_raw(legs{i}(2), :)')
%     plot(time, gcomp_imu(legs{i}(2), :)')
    
    subplot(3,3,7)
%     plot(time, joint_effort(legs{i}(7), :)')
    xlim(xlims)
    xlabel('Time [s]')
    ylabel('Torque [Nm]')
    title([leg_names{i}, ' 7'])
    hold on
    grid on
    plot(time, effort_reference(legs{i}(7), :)'- tau_m_filt(legs{i}(7), :)', 'LineWidth', 3)
%     plot(time, effort_reference(legs{i}(7), :)')
%     plot(time, gcomp_raw(legs{i}(3), :)')
%     plot(time, gcomp_imu(legs{i}(3), :)')
    
%     fig = figure;
%     fig.Units = 'centimeters';
%     fig.Position = [10 10 25 17];
%     fig.Color = [1 1 1];
%     
%     plot(time, tau_m_filt(legs{i}(1:4),:)' - gcomp_imu(legs{i}(1:4),:)', 'LineWidth', 2)
%     xlim(xlims)
%     xlabel('Time [s]')
%     ylabel('Torque [Nm]')
%     title([leg_names{i}, ' Torque offset']);
%     hold on
%     grid on
%     legend({' Hip pitch', ' Hip pitch', ' Hip yaw', ' Knee pitch'}, 'Location', 'best')
    
    savefig(['/tmp/torque_offsets_'  leg_names{i} '.fig']);
    print(['/tmp/torque_offsets_'  leg_names{i} '.eps'], '-depsc');
    

end