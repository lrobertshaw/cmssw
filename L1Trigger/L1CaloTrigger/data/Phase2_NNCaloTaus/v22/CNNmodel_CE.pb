
S
TowerClusterImagePlaceholder*
dtype0*$
shape:���������	
P
AssociatedCl3dFeaturesPlaceholder*
dtype0*
shape:���������
N
TowerClusterPositionPlaceholder*
dtype0*
shape:���������
�
;TauMinator_CE_conv/CNNlayer1/Conv2D/ReadVariableOp/resourceConst*
dtype0*�
value�B�"�I��>���J8�>�N����>Q>w.�>9f�>3�>�I���:�=��>9Y�Ԭ ��$x>�A?�C�YZa>�՗>�7J>j9�=Su��TZ0?=>/S	><����>hUt>�Sܾ'-��l�C>�6~=
�
2TauMinator_CE_conv/CNNlayer1/Conv2D/ReadVariableOpIdentity;TauMinator_CE_conv/CNNlayer1/Conv2D/ReadVariableOp/resource*
T0
�
#TauMinator_CE_conv/CNNlayer1/Conv2DConv2DTowerClusterImage2TauMinator_CE_conv/CNNlayer1/Conv2D/ReadVariableOp*
T0*
data_formatNHWC*
	dilations
*
explicit_paddings
 *
paddingVALID*
strides
*
use_cudnn_on_gpu(
t
7TauMinator_CE_conv/BN_CNNlayer1/ReadVariableOp/resourceConst*
dtype0*%
valueB"�R?l�?���?�|x?
|
.TauMinator_CE_conv/BN_CNNlayer1/ReadVariableOpIdentity7TauMinator_CE_conv/BN_CNNlayer1/ReadVariableOp/resource*
T0
v
9TauMinator_CE_conv/BN_CNNlayer1/ReadVariableOp_1/resourceConst*
dtype0*%
valueB"�7�>��o�*���>
�
0TauMinator_CE_conv/BN_CNNlayer1/ReadVariableOp_1Identity9TauMinator_CE_conv/BN_CNNlayer1/ReadVariableOp_1/resource*
T0
�
HTauMinator_CE_conv/BN_CNNlayer1/FusedBatchNormV3/ReadVariableOp/resourceConst*
dtype0*%
valueB" (�>*ch�8��?o'J?
�
?TauMinator_CE_conv/BN_CNNlayer1/FusedBatchNormV3/ReadVariableOpIdentityHTauMinator_CE_conv/BN_CNNlayer1/FusedBatchNormV3/ReadVariableOp/resource*
T0
�
JTauMinator_CE_conv/BN_CNNlayer1/FusedBatchNormV3/ReadVariableOp_1/resourceConst*
dtype0*%
valueB"!@Y��@�X�@Rb@
�
ATauMinator_CE_conv/BN_CNNlayer1/FusedBatchNormV3/ReadVariableOp_1IdentityJTauMinator_CE_conv/BN_CNNlayer1/FusedBatchNormV3/ReadVariableOp_1/resource*
T0
�
0TauMinator_CE_conv/BN_CNNlayer1/FusedBatchNormV3FusedBatchNormV3#TauMinator_CE_conv/CNNlayer1/Conv2D.TauMinator_CE_conv/BN_CNNlayer1/ReadVariableOp0TauMinator_CE_conv/BN_CNNlayer1/ReadVariableOp_1?TauMinator_CE_conv/BN_CNNlayer1/FusedBatchNormV3/ReadVariableOpATauMinator_CE_conv/BN_CNNlayer1/FusedBatchNormV3/ReadVariableOp_1*
T0*
U0*
data_formatNHWC*
epsilon%o�:*
exponential_avg_factor%  �?*
is_training( 
i
&TauMinator_CE_conv/RELU_CNNlayer1/ReluRelu0TauMinator_CE_conv/BN_CNNlayer1/FusedBatchNormV3*
T0
�
'TauMinator_CE_conv/MP_CNNlayer1/MaxPoolMaxPool&TauMinator_CE_conv/RELU_CNNlayer1/Relu*
T0*
data_formatNHWC*
explicit_paddings
 *
ksize
*
paddingVALID*
strides

�
;TauMinator_CE_conv/CNNlayer2/Conv2D/ReadVariableOp/resourceConst*
dtype0*�
value�B�"����6��s=q�I�>���=�8�>1�4��*k��߁=4+�}�����h�<�9�>���<���<��Z�;��J�=kھ|;��q�P�#�>~�<3[ݽ��]�Й�=�
�f����!L�	윽�jٽ�@��Q��Ȟ���j>�����(�=,����=y>0�V>��=�G;��Y=P���̽��=B����z>���=m�[-*��c�=87�<>ަ���"=z!}�f@�>�@}�'Q5<Ƣ���x�0�D>�
��Z�Ծ�����h=N�=7n�_I�>k���+W�B��=Щ�o�	>b��>��z>����<>����-eݽ��>Ջ%���¾*�>#�$>��Ⱦy0��l��5�@=q3����>Y�=R���ȅ����������"�<�N>�l>�R�JS�� w���q>_����=��?� �0��'?����OC>�Ն��F>�T�Lͽ��;�?���>���p���1�������>�yL=u�>�/:
�
2TauMinator_CE_conv/CNNlayer2/Conv2D/ReadVariableOpIdentity;TauMinator_CE_conv/CNNlayer2/Conv2D/ReadVariableOp/resource*
T0
�
#TauMinator_CE_conv/CNNlayer2/Conv2DConv2D'TauMinator_CE_conv/MP_CNNlayer1/MaxPool2TauMinator_CE_conv/CNNlayer2/Conv2D/ReadVariableOp*
T0*
data_formatNHWC*
	dilations
*
explicit_paddings
 *
paddingVALID*
strides
*
use_cudnn_on_gpu(
�
7TauMinator_CE_conv/BN_CNNlayer2/ReadVariableOp/resourceConst*
dtype0*5
value,B*" ��?��?k��?�`?�~�?��?�`?䝋?
|
.TauMinator_CE_conv/BN_CNNlayer2/ReadVariableOpIdentity7TauMinator_CE_conv/BN_CNNlayer2/ReadVariableOp/resource*
T0
�
9TauMinator_CE_conv/BN_CNNlayer2/ReadVariableOp_1/resourceConst*
dtype0*5
value,B*" �i�d�=��>'�>f:�>(�>�+�>��>
�
0TauMinator_CE_conv/BN_CNNlayer2/ReadVariableOp_1Identity9TauMinator_CE_conv/BN_CNNlayer2/ReadVariableOp_1/resource*
T0
�
HTauMinator_CE_conv/BN_CNNlayer2/FusedBatchNormV3/ReadVariableOp/resourceConst*
dtype0*5
value,B*" Oc������fֿ����_��o�`�Ko�>�:��
�
?TauMinator_CE_conv/BN_CNNlayer2/FusedBatchNormV3/ReadVariableOpIdentityHTauMinator_CE_conv/BN_CNNlayer2/FusedBatchNormV3/ReadVariableOp/resource*
T0
�
JTauMinator_CE_conv/BN_CNNlayer2/FusedBatchNormV3/ReadVariableOp_1/resourceConst*
dtype0*5
value,B*" ¾�?��9@���?���?�)�@}]�?�G?L�?
�
ATauMinator_CE_conv/BN_CNNlayer2/FusedBatchNormV3/ReadVariableOp_1IdentityJTauMinator_CE_conv/BN_CNNlayer2/FusedBatchNormV3/ReadVariableOp_1/resource*
T0
�
0TauMinator_CE_conv/BN_CNNlayer2/FusedBatchNormV3FusedBatchNormV3#TauMinator_CE_conv/CNNlayer2/Conv2D.TauMinator_CE_conv/BN_CNNlayer2/ReadVariableOp0TauMinator_CE_conv/BN_CNNlayer2/ReadVariableOp_1?TauMinator_CE_conv/BN_CNNlayer2/FusedBatchNormV3/ReadVariableOpATauMinator_CE_conv/BN_CNNlayer2/FusedBatchNormV3/ReadVariableOp_1*
T0*
U0*
data_formatNHWC*
epsilon%o�:*
exponential_avg_factor%  �?*
is_training( 
i
&TauMinator_CE_conv/RELU_CNNlayer2/ReluRelu0TauMinator_CE_conv/BN_CNNlayer2/FusedBatchNormV3*
T0
Z
%TauMinator_CE_conv/CNNflattened/ConstConst*
dtype0*
valueB"����   
�
'TauMinator_CE_conv/CNNflattened/ReshapeReshape&TauMinator_CE_conv/RELU_CNNlayer2/Relu%TauMinator_CE_conv/CNNflattened/Const*
T0*
Tshape0
R
(TauMinator_CE_conv/middleMan/concat/axisConst*
dtype0*
value	B :
�
#TauMinator_CE_conv/middleMan/concatConcatV2'TauMinator_CE_conv/CNNflattened/ReshapeTowerClusterPositionAssociatedCl3dFeatures(TauMinator_CE_conv/middleMan/concat/axis*
N*
T0*

Tidx0
�
NoOpNoOp@^TauMinator_CE_conv/BN_CNNlayer1/FusedBatchNormV3/ReadVariableOpB^TauMinator_CE_conv/BN_CNNlayer1/FusedBatchNormV3/ReadVariableOp_1/^TauMinator_CE_conv/BN_CNNlayer1/ReadVariableOp1^TauMinator_CE_conv/BN_CNNlayer1/ReadVariableOp_1@^TauMinator_CE_conv/BN_CNNlayer2/FusedBatchNormV3/ReadVariableOpB^TauMinator_CE_conv/BN_CNNlayer2/FusedBatchNormV3/ReadVariableOp_1/^TauMinator_CE_conv/BN_CNNlayer2/ReadVariableOp1^TauMinator_CE_conv/BN_CNNlayer2/ReadVariableOp_13^TauMinator_CE_conv/CNNlayer1/Conv2D/ReadVariableOp3^TauMinator_CE_conv/CNNlayer2/Conv2D/ReadVariableOp*"
_acd_function_control_output(
I
IdentityIdentity#TauMinator_CE_conv/middleMan/concat^NoOp*
T0"�