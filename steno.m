disp('Welcome to the Image Encryption and Steganography Using Contourlet Transform and Chaotic Maps demo!');

disp('Here is the process of the first part of program');
disp('Encrypting and hiding image:')
disp('1. Plain secret image (Lena) is read and encrypt using Contourlet Transform and Chaotic Maps.')
disp('2. Embed (hide) encrypted image into Cover image (Sydney.png).')
disp(' ');

%input( 'Press Enter key to continue...' ) ;
disp( ' ' );

% Constants for chaotic function
initial_value = 0.815647;
r_param = 3.987;

% Test images
secret_img_name = 'lena.png';
cover_img_name  = 'sydney_1024.png';

% images converted to grayscale
secret_image = imread(secret_img_name);
if size(secret_image, 3) == 3
    secret_image = rgb2gray(secret_image);  
end   

% Encrypt the image
encrypted_image = ContourletEncDec.chaotic_image_encryption(secret_image);

% Display Original-Encryted-Decryted images
disp('Display Original-Encryted-Decryted images');
figure;
subplot(1, 2, 1); imshow(secret_image);             title('Original Image');
subplot(1, 2, 2); imshow(uint8(encrypted_image));   title('Encrypted Image');

%input( 'Press Enter key to continue...' ) ;
disp( ' ' );

cover_image = imread(cover_img_name);
if size(cover_image, 3) == 3
    cover_image = rgb2gray(cover_image);  
end   

% embedding process
stego_image= embed_image(cover_image, encrypted_image, initial_value, r_param);
%imwrite(stego_image, 'stego_image.png');

% Cover image + Secret image -> Embedded image
disp('Display Cover image + Secret image -> Embedded image');
figure;
subplot(1,3,1);imshow(cover_image);title('Cover Image');
subplot(1,3,2);imshow(encrypted_image);title('Secret Image');
subplot(1,3,3);imshow(stego_image);title('Embedded Image');

%input( 'Press Enter key to continue...' ) ;
disp( ' ' );

disp('Here is the process of the second part of the program');
disp('Extracting secret image and decrypt secret image:')
disp('1. Extract the secret image from Cover image (Sydney)')
disp('2. Secret image is decrypted using Contourlet Transform and Chaotic Maps with inverse orders.')
disp(' ');

extracted_image = extract_image(stego_image, size(encrypted_image), initial_value, r_param);
%imwrite(uint8(extracted_image), 'extracted_image.png');


% Secret image vs extracted image
disp('Display Secret image vs extracted image');
figure;
subplot(1,2,1);imshow(encrypted_image);title('Secret Image');
subplot(1,2,2);imshow(uint8(extracted_image));title('Extracted Secret Image');

%input( 'Press Enter key to continue...' ) ;

% Decrypt the extracted image
decrypted_image = ContourletEncDec.chaotic_image_decryption(extracted_image);

% Display Original-Decrypted images
disp('Display Original-Decrypted images');
figure;
subplot(1, 2, 1); imshow(uint8(extracted_image));   title('Extracted Image');
subplot(1, 2, 2); imshow(uint8(decrypted_image));   title('Decrypted Image');

%input( 'Press Enter key to continue...' ) ;


disp('Display STEGOANALYSIS');

% 1. Histogram Analysis
figure;
subplot(1, 2, 1); imhist(uint8(cover_image)); title('Cover Image Histogram');
subplot(1, 2, 2); imhist(uint8(stego_image)); title('Stego Image Histogram');

% 2. Calculate PSNR 
psnr_val = psnr(stego_image,cover_image);
fprintf('PSNR between cover and stego image: %.2f dB\n', psnr_val);

% 3. Calculate SSIM
ssim_val = ssim(stego_image, cover_image);
fprintf('SSIM between cover and stego image: %.2f dB\n', ssim_val);

% 4. Chi-Square Analysis 
chi_square_analysis(uint8(cover_image),uint8(stego_image));

% 5. NCC
ncc_analysis(cover_image,stego_image);

% Embeds a secret image into a cover image using the Least Significant Bit (LSB) method.
%
% Args:
%   cover_image: The cover image (uint8).
%   secret_image: The secret image (uint8).
%   chaotic params : Chaotic sequence parameters (initial_value and r_param)
%
% Returns:
%   stego_image: The stego image (uint8) containing the hidden secret image.
function embedded_image = embed_image(cover_image, secret_image, initial_value, r_param)

    secret_binary = dec2bin(secret_image(:), 8);
    secret_bits = reshape(secret_binary', 1, []);
    num_bits_to_embed = length(secret_bits);

    % Pre-calculate chaotic sequence indices
    num_elements = numel(cover_image);
    % chaotic sequence içinde tekilliği sağlamak için daha fazla sequence
    % oluşturuluyor
    chaos_sequence = generate_chaotic_sequence(num_elements*3, initial_value, r_param);
    chaos_sequence_indices = unique(mod(floor(chaos_sequence * 10^15), num_elements) + 1);

    embedded_image = cover_image;
    bit_index = 1;    
    for k = 1:num_bits_to_embed
        index = chaos_sequence_indices(mod(k-1,length(chaos_sequence_indices))+1); % Use mod to wrap around indices
        [row, col, channel] = ind2sub(size(cover_image), index);
        
        cover_pixel = uint8(embedded_image(row, col, channel));
        secret_bit = str2double(secret_bits(k));
        embedded_image(row, col, channel) = bitset(cover_pixel, 1, secret_bit);
    end
end

% Extacts a secret image from a steno image using the Least Significant Bit (LSB) method.
%
% Args:
%   cover_image: Steno image (uint8).
%   secret_image_size: The secret image size 
%   chaotic params : Chaotic sequence parameters (initial_value and r_param)
%
% Returns:
%   stego_image: The stego image (uint8) containing the hidden secret image.
function extracted_image = extract_image(embedded_image, secret_size, initial_value, r_param)

    num_pixels = prod(secret_size);
    num_bits_to_extract = num_pixels * 8;
    
    % Pre-calculate chaotic sequence indices
    num_elements = numel(embedded_image);
    % chaotic sequence içinde tekilliği sağlamak için daha fazla sequence
    % oluşturuluyor
    chaos_sequence = generate_chaotic_sequence(num_elements*3, initial_value, r_param);
    chaos_sequence_indices = unique(mod(floor(chaos_sequence * 10^15), num_elements) + 1);

    extracted_bits = '';
    
    for k = 1:num_bits_to_extract
        index = chaos_sequence_indices(mod(k-1,length(chaos_sequence_indices))+1); % Use mod to wrap around indices
        [row, col, channel] = ind2sub(size(embedded_image), index);

        extracted_bits = [extracted_bits num2str(bitget(uint8(embedded_image(row, col, channel)), 1))];
    end

    extracted_binary = reshape(extracted_bits, 8, []);
    extracted_decimal = bin2dec(extracted_binary');
    extracted_image = reshape(extracted_decimal, secret_size);
end


% Chaotic Sequence Generation (Logistic Map)
% Logistic map: x_(n+1) = r * x_n * (1 - x_n)
function chaos_sequence = generate_chaotic_sequence(length, initial_value, r_param)
    chaos_sequence = zeros(1, length);
    chaos_sequence(1) = initial_value;
    %r_param = 3.987;
    for i = 2:length
        chaos_sequence(i) = r_param * chaos_sequence(i-1) * (1 - chaos_sequence(i-1));
    end
end


function chi_square_analysis(cover_image, stego_image)

    % Histogram
    cover_hist = imhist(cover_image);
    stego_hist = imhist(stego_image);

    % Chi-square test a
    [chi2_stat, p_value] = chi2gof(stego_hist, 'Expected', cover_hist);

    % Results
    fprintf('Chi-Kare İstatistiği: %f\n', chi2_stat);
    fprintf('P-Değeri: %f\n', p_value);

end

function ncc_analysis(cover_image, stego_image)

    % NCC hesapla
    ncc_value = normxcorr2(double(stego_image), double(cover_image));

    % NCC matrisinin ortasındaki değeri al (tam örtüşme noktası)
    [max_ncc, ~] = max(abs(ncc_value(:)));

    % Görüntüleri ve NCC sonucunu göster
    figure;    
    surf(ncc_value), shading flat;    
    zlabel('NCC Değeri');

    fprintf('Maksimum NCC Değeri: %f\n', max_ncc);
end