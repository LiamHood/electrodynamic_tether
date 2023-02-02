function [I, debug_values] = controlled_current(limit_libration, current_value, roll, pitch, droll, dpitch, m, L)
%     lim_e = limit_libration;    
%     I = current_value;
%     pe_roll = roll;
%     ke_roll = droll^3/abs(droll)*6e4;
%     e_roll = pe_roll+ke_roll;
%     pe_pitch = pitch;
%     ke_pitch = dpitch^3/abs(dpitch)*6e4;
%     e_pitch = pe_pitch+ke_pitch;
%     if abs(e_roll) > lim_e
%         I = 0;
%     elseif abs(e_pitch) > lim_e
%         I = 0;
%     end

    I = current_value;
    lim_e = L*(1-cos(deg2rad(limit_libration)));
    pe_roll = L*(1-cos(roll));
    ke_roll = (1/2)*(L*droll)^2;
    if roll < 0                   
        pe_roll = -pe_roll;
    end
    if droll < 0
        ke_roll = -ke_roll;
    end
    e_roll = pe_roll + ke_roll;
    pe_pitch = L*(1-cos(pitch));
    ke_pitch = (1/2)*(L*dpitch)^2;
%             ke_t = (1/2)*Iy*(dtheta)^2;

    if pitch < 0
        pe_pitch = -pe_pitch;
    end
    if dpitch < 0
        ke_pitch = -ke_pitch;
    end
    e_pitch = pe_pitch + ke_pitch;
    if abs(e_roll) > lim_e
        I = 0;
    elseif abs(e_pitch) > lim_e
        I = 0;
    end


    debug_values = [lim_e, pe_roll, ke_roll, e_roll, pe_pitch, ke_pitch, e_pitch];
end