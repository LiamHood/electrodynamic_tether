function C_attitude = tether_attitude(in_plane, out_plane)
    u1 = [cos(out_plane)*cos(in_plane); -sin(out_plane); cos(out_plane)*sin(in_plane)];
    u2 = [-sin(in_plane); 0; cos(in_plane)];
    u3 = [-sin(out_plane)*cos(in_plane); -cos(out_plane); -sin(out_plane)*sin(in_plane)];
    C_attitude = [u1, u2, u3];
end